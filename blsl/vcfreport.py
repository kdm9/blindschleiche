#!/usr/bin/env python
import plotly.express as px
from cyvcf2 import VCF
import pandas as pd
import numpy as np
from tqdm.auto import tqdm
import physt
from physt import h1 as histogram
import physt.plotting
physt.plotting.set_default_backend("plotly")
import sklearn.impute
import sklearn.decomposition


from sys import stdin, stdout, stderr
from subprocess import Popen, PIPE, DEVNULL
import argparse
import math
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
import random


def parallel_regions(vcf, chunk=1000000):
    V = VCF(vcf)
    for cname, clen in zip(V.seqnames, V.seqlens):
        for start in range(0, clen, chunk):
            s = start+1
            e = min(start+chunk, clen)
            yield f"{cname}:{s}-{e}"



def variant2dict(v, fields=None, genotypes=False):
    dat = {"CHROM": v.CHROM, "POS": v.POS, "REF": v.REF, "ALT": v.ALT, "QUAL": v.QUAL}
    dat["call_rate"] = v.call_rate
    for key, val  in v.INFO:
        if isinstance(val, str) and "," in val:
            val = val.split(',')
        if isinstance(val, tuple) or isinstance(val, list):
            for i, alt in enumerate(v.ALT):
                dat[f"INFO_{key}_{i}"] = val[i]
        else:
            dat[f"INFO_{key}"] = val
    if fields:
        dat = {K:V for K, V in dat.items() if K in fields}
    if genotypes:
        dat["genotypes"] = [v[0] + v[1] if v[0] != -1 and v[1] != -1 else float('nan') for v in v.genotypes ]
    return dat


def one_chunk_stats(vcf, chunk, fill=True, fields=None, min_maf=0.01, min_call=0.7, subsample=0.01):
    cmd=f"bcftools view -r {chunk} {vcf} -Ou"
    if fill:
        cmd = f"{cmd} | bcftools +fill-tags - -Ou -- -d -t all,F_MISSING"
    variants = []
    res = []
    with Popen(cmd, shell=True, stdout=PIPE, stderr=DEVNULL) as proc:
        vcf = VCF(proc.stdout)
        for v in vcf:
            if v.call_rate >= min_call and len(v.ALT) == 1 and v.INFO["MAF"] >= min_maf and random.random() <= subsample:
                variants.append(variant2dict(v, genotypes=True))
            res.append(variant2dict(v, fields=["INFO_MAF", "call_rate", "INFO_HWE", "INFO_ExcHet", "INFO_AC", "QUAL"]))
        del vcf
    if not res:
        return {
            "chunk": chunk,
            "n_snps": len(res),
        }
    df = pd.DataFrame.from_records(res)
    try:
        return {
            "chunk": chunk,
            "n_snps": len(res),
            "maf":  histogram(df.INFO_MAF, "fixed_width", bin_width=0.01, adaptive=True),
            "call_rate":  histogram(df.call_rate, "fixed_width", bin_width=0.01,  adaptive=True),
            "hwe":  histogram(df.INFO_HWE, "fixed_width", bin_width=0.01, adaptive=True),
            "exc_het": histogram(df.INFO_ExcHet, "fixed_width", bin_width=0.01, adaptive=True),
            "ac": histogram(df.INFO_AC, "fixed_width", bin_width=1, adaptive=True),
            "qual": histogram(df.QUAL, "fixed_width", bin_width=1, adaptive=True),
            "subset_variants": variants
        }
    except Exception as exc:
        print(exc.__class__.__name__, str(exc), df, res)



def update_result(globalres, res, hists = ["maf", "call_rate", "hwe", "exc_het", "ac", "qual"], sums=["n_snps"], lists=["subset_variants",]):
    for listkey in lists:
        if listkey not in res:
            continue
        if listkey not in globalres:
            globalres[listkey] = []
        globalres[listkey].extend(res[listkey])
    for sumkey in sums:
        if sumkey not in res:
            continue
        if sumkey not in globalres:
            globalres[sumkey] = 0
        globalres[sumkey] += res[sumkey]
    for histkey in hists:
        if histkey not in res:
            continue
        if histkey not in globalres:
            globalres[histkey] = res[histkey]
            continue
        globalres[histkey] += res[histkey]
    return globalres


def chunkwise_bcftools_stats(vcf, threads=8):
    regions = list(parallel_regions(vcf))
    global_res = {
        "samples": VCF(vcf).samples
    }
    with tqdm(total=len(regions), unit="chunk") as pbar:
        for i in range(0, len(regions), 1000):
            to=min(len(regions), i+1000)
            with ProcessPoolExecutor(threads) as exc:
                jobs = (exc.submit(one_chunk_stats, vcf, region) for region in regions[i:to])
                for job in as_completed(jobs):
                    pbar.update(1)
                    update_result(global_res, job.result())
    return global_res


def genotype_pca(res):
    genomat = np.vstack([x["genotypes"] for x in res])
    imp = sklearn.impute.SimpleImputer()
    genomat = imp.fit_transform(genomat)
    pc = sklearn.decomposition.PCA()
    gpc = pc.fit_transform(genomat.T)
    return genomat, gpc, pc.explained_variance_ratio_


def generate_report(vcf, threads):
    gres = chunkwise_bcftools_stats(vcf, threads=threads)
    
    fig_maf = gres["maf"].plot()
    fig_maf.update_xaxes(title_text="MAF")
    fig_maf.update_yaxes(title_text="# SNPS")
    fig_maf.update_layout(title_text="MAF Spectrum")
    MAF_CODE=fig_maf.to_html(full_html=False)

    fig_cr = gres["call_rate"].plot()
    fig_cr.update_xaxes(title_text="Call Rate")
    fig_cr.update_yaxes(title_text="# SNPS")
    fig_cr.update_layout(title_text="Call Rate Spectrum")
    CALLRATE_CODE = fig_cr.to_html(full_html=False)

    fig_hwe = gres["hwe"].plot()
    fig_hwe.update_xaxes(title_text="p(not in HWE)")
    fig_hwe.update_yaxes(title_text="# SNPS")
    fig_hwe.update_layout(title_text="HWE Test")
    HWE_CODE = fig_hwe.to_html(full_html=False)
    
    fig_xh = gres["exc_het"].plot()
    fig_xh.update_xaxes(title_text="p(Excess Heterozygotes)")
    fig_xh.update_yaxes(title_text="# SNPS")
    fig_xh.update_layout(title_text="ExHet Test")
    EXHET_CODE = fig_xh.to_html(full_html=False)

    fig_ac = gres["ac"].plot()
    fig_ac.update_xaxes(title_text="Alternate Allele Count")
    fig_ac.update_yaxes(title_text="# SNPS")
    fig_ac.update_layout(title_text="Count of Alt Alleles")
    AC_CODE = fig_ac.to_html(full_html=False)
    
    fig_qual = gres["qual"].plot()
    fig_qual.update_xaxes(title_text="Variant Quality")
    fig_qual.update_yaxes(title_text="# SNPS")
    fig_qual.update_layout(title_text="Variant Quality")
    QUAL_CODE = fig_qual.to_html(full_html=False)

    
    x, pc, varex = genotype_pca(gres["subset_variants"])
    pc_fig = px.scatter(
        x=pc[:,0],
        y=pc[:,1],
        labels={"xy"[i]: f"PC {i+1} ({varex[i]*100:.1f}%)" for i in range(2)},
        hover_name=gres["samples"],
        title="PCA on imputed subset of high-quality SNPs",
    )
    PCA_CODE = pc_fig.to_html(full_html=False)
    
    html = f"""
    <html>
    <head>
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <meta http-equiv="X-UA-Compatible" content="ie=edge">
    </head>
    <body>
    <h1>Statistics</h1>
    <h2>Minor Allele Frequency</h2>
    {MAF_CODE}
    <h2>Call Rate</h2>
    {CALLRATE_CODE}
    <h2>Hardy-Weinberg Equilibirum Tests</h2>
    {HWE_CODE}
    <h2>Excess Heterozygosity Tests</h2>
    {EXHET_CODE}
    <h2>Allele Counts</h2>
    {AC_CODE}
    <h2>Variant Quality</h2>
    {QUAL_CODE}
    <h1>PCA</h1>
    {PCA_CODE}
    </body>
    </html>
    """
    return html


def main(argv=None):
    """vcfreport: Prepare a basic html report about a VCF file"""
    ap = argparse.ArgumentParser("blsl vcfreport")
    ap.add_argument("--output", "-o", type=argparse.FileType("w"), required=True,
                    help="Output html file")
    ap.add_argument("--threads", "-j", type=int, default=2,
                    help="Parallel threads")
    ap.add_argument("vcf",
                    help="VCF input file (must be indexed)")
    args=ap.parse_args(argv)

    html = generate_report(args.vcf, threads=args.threads)
    args.output.write(html)
    args.output.flush()
