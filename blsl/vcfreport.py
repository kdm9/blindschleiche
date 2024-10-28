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
from dataclasses import dataclass


@dataclass
class Region:
    contig: str
    start: int
    end: int
    gpos: int

    def __str__(self):
        return f"{self.contig}:{self.start}-{self.end}"

    @property
    def len(self):
        return self.end - self.start + 1


def parallel_regions(vcf, chunk=1000000):
    V = VCF(vcf)
    gclen = 0
    for cname, clen in zip(V.seqnames, V.seqlens):
        for start in range(0, clen, chunk):
            s = start+1
            e = min(start+chunk, clen)
            gpos = gclen + ( s + e )/2
            yield Region(cname, s, e, gpos)
        gclen += clen


def variant2dict(v, fields=None):
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
    if "FORMAT_GT" in fields:
        dat["FORMAT_GT"] = [v[0] + v[1] if v[0] != -1 and v[1] != -1 else float('nan') for v in v.genotypes ]
    if "FORMAT_DP" in fields:
        dat["FORMAT_DP"] = v.format('DP')
    return dat


def one_chunk_stats(vcf, chunk, fill=True, fields=None, min_maf=0.01, min_call=0.7, subsample=0.01):
    cmd=f"bcftools view -r {str(chunk)} {vcf} -Ou"
    if fill:
        cmd = f"{cmd} | bcftools +fill-tags - -Ou -- -d -t all,F_MISSING"
    variants = []
    res = []
    with Popen(cmd, shell=True, stdout=PIPE, stderr=DEVNULL) as proc:
        vcf = VCF(proc.stdout)
        for v in vcf:
            if v.call_rate >= min_call and len(v.ALT) == 1 and v.INFO["MAF"] >= min_maf and random.random() <= subsample:
                variants.append(variant2dict(v, fields=["FORMAT_GT", "FORMAT_DP"]))
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
            "ac": histogram(df.INFO_AC, "fixed_width", bin_width=10, adaptive=True),
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
    if "snps_per_window" not in globalres:
        globalres["snps_per_window"] = {"pos": list(), "n_snps": list(), "chrom": list()}
    globalres["snps_per_window"]["pos"].append(res["chunk"].gpos)
    globalres["snps_per_window"]["chrom"].append(res["chunk"].contig)
    globalres["snps_per_window"]["n_snps"].append(res["n_snps"]/res["chunk"].len)
    return globalres


def chunkwise_bcftools_stats(vcf, threads=8, chunksize=1_000_000):
    regions = list(parallel_regions(vcf, chunk=chunksize))
    v=VCF(vcf)
    global_res = {
        "samples": v.samples,
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

def sample_missingness(variants):
    genomat = np.vstack([x["FORMAT_GT"] for x in variants])
    missing = np.sum(np.isnan(genomat), axis=0) / np.shape(genomat)[0]
    #return {s:m for s, m in zip(res["samples"], missing)}
    return histogram(missing, 30)

def sample_depth(variants):
    dpmat = np.hstack([x["FORMAT_DP"] for x in variants])
    meandepth = np.nansum(dpmat, axis=1) / np.shape(dpmat)[1]
    #return {s:m for s, m in zip(samples, meandepth)}
    return histogram(meandepth, 30)

def genotype_pca(res):
    genomat = np.vstack([x["FORMAT_GT"] for x in res])
    imp = sklearn.impute.SimpleImputer()
    genomat = imp.fit_transform(genomat)
    pc = sklearn.decomposition.PCA()
    gpc = pc.fit_transform(genomat.T)
    return genomat, gpc, pc.explained_variance_ratio_


def generate_report(vcf, threads, chunksize=1_000_000):
    gres = chunkwise_bcftools_stats(vcf, threads=threads, chunksize=chunksize)
    figwidth=1000
    figheight=750
    
    fig_smis = sample_missingness(gres["subset_variants"]).plot()
    fig_smis.update_xaxes(title_text="Missing Rate (sample)")
    fig_smis.update_yaxes(title_text="# Samples")
    fig_smis.update_layout(title_text="Sample Missing Rate", width=figwidth, height=figheight)
    SMIS_CODE=fig_smis.to_html(full_html=False)

    fig_sdp = sample_depth(gres["subset_variants"]).plot()
    fig_sdp.update_xaxes(title_text="Mean Depth (sample)")
    fig_sdp.update_yaxes(title_text="# Samples")
    fig_sdp.update_layout(title_text="Sample Mean Depths", width=figwidth, height=figheight)
    SDP_CODE=fig_sdp.to_html(full_html=False)

    fig_maf = gres["maf"].plot()
    fig_maf.update_xaxes(title_text="MAF")
    fig_maf.update_yaxes(title_text="# SNPS")
    fig_maf.update_layout(title_text="MAF Spectrum", width=figwidth, height=figheight)
    MAF_CODE=fig_maf.to_html(full_html=False)

    fig_cr = gres["call_rate"].plot()
    fig_cr.update_xaxes(title_text="Call Rate")
    fig_cr.update_yaxes(title_text="# SNPS")
    fig_cr.update_layout(title_text="Call Rate Spectrum", width=figwidth, height=figheight)
    CALLRATE_CODE = fig_cr.to_html(full_html=False)

    fig_hwe = gres["hwe"].plot()
    fig_hwe.update_xaxes(title_text="p(not in HWE)")
    fig_hwe.update_yaxes(title_text="# SNPS")
    fig_hwe.update_layout(title_text="HWE Test", width=figwidth, height=figheight)
    HWE_CODE = fig_hwe.to_html(full_html=False)
    
    fig_xh = gres["exc_het"].plot()
    fig_xh.update_xaxes(title_text="p(Excess Heterozygotes)")
    fig_xh.update_yaxes(title_text="# SNPS")
    fig_xh.update_layout(title_text="ExHet Test", width=figwidth, height=figheight)
    EXHET_CODE = fig_xh.to_html(full_html=False)

    fig_ac = gres["ac"].plot()
    fig_ac.update_xaxes(title_text="Alternate Allele Count")
    fig_ac.update_yaxes(title_text="# SNPS")
    fig_ac.update_layout(title_text="Count of Alt Alleles", width=figwidth, height=figheight)
    AC_CODE = fig_ac.to_html(full_html=False)
    
    fig_qual = gres["qual"].plot()
    fig_qual.update_xaxes(title_text="Variant Quality")
    fig_qual.update_yaxes(title_text="# SNPS")
    fig_qual.update_layout(title_text="Variant Quality", width=figwidth, height=figheight)
    QUAL_CODE = fig_qual.to_html(full_html=False)

    
    x, pc, varex = genotype_pca(gres["subset_variants"])
    axtitle={"xy"[i]: f"PC {i+1} ({varex[i]*100:.1f}%)" for i in range(2)}
    pc_fig = px.scatter(
        x=pc[:,0],
        y=pc[:,1],
        hover_name=gres["samples"],
        width=figwidth,
        height=figheight,
    )
    pc_fig.update_xaxes(title_text=axtitle["x"])
    pc_fig.update_yaxes(title_text=axtitle["y"])
    pc_fig.update_layout(title_text="PCA on imputed subset of high-quality SNPs", width=figwidth, height=figheight)
    PCA_CODE = pc_fig.to_html(full_html=False)
    
    sog_data = pd.DataFrame(gres["snps_per_window"])
    sog_data = sog_data.sort_values(by="pos")
    sog_fig = px.line(
        sog_data,
        x="pos",
        y="n_snps",
        color="chrom",
        title="SNPs per Genome Window",
        width=figwidth,
        height=figheight,
    )
    sog_fig.update_xaxes(title_text="Genome Position")
    sog_fig.update_yaxes(title_text="SNPs per base")
    sog_fig.update_layout(title_text="SNPs per Genome Window", width=figwidth, height=figheight)
    SNPOVERGENOME_CODE = sog_fig.to_html(full_html=False)

    html = f"""
    <html>
    <head>
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <meta http-equiv="X-UA-Compatible" content="ie=edge">
        <style>
            body {{
              max-width: 1000px;
            }}
        </style>
    </head>
    <body>
        <h1>VCF Statistics</h1>
        <table>
            <tr><td>File</td> <td align="right">{vcf}</td> </tr>
            <tr><td>Total # SNPs</td> <td align="right">{gres['n_snps']:,}</td> </tr>
            <tr><td># Samples</td> <td align="right">{len(gres['samples']):,}</td> </tr>
        </table>

        <h1>Sample-level Stats</h1>
        <p>These statstics are based on a random 1% of all SNPs for efficiency's sake</p>

        <h2>Sample Missing Rate</h1>
        {SMIS_CODE}

        <h2>Sample Mean Depth</h1>
        {SDP_CODE}

        <h2>Sample PCA</h1>
        {PCA_CODE}

        <h1> Variant-level stats</h1>
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

        <h2>SNPs per Genome Window</h1>
        {SNPOVERGENOME_CODE}

    </body>
    </html>
    """
    return html


def main(argv=None):
    """vcfreport: Prepare a basic html report about a VCF file"""
    ap = argparse.ArgumentParser("blsl vcfreport")
    ap.add_argument("--output", "-o", type=argparse.FileType("w"), required=True,
                    help="Output html file")
    ap.add_argument("--threads", "-t", type=int, default=2,
                    help="Parallel threads")
    ap.add_argument("--chunksize", "-c", type=int, default=1_000_000,
                    help="Chunks size for parallelism")
    ap.add_argument("vcf",
                    help="VCF input file (must be indexed)")
    args=ap.parse_args(argv)

    html = generate_report(args.vcf, threads=args.threads, chunksize=args.chunksize)
    args.output.write(html)
    args.output.flush()

if __name__ == "__main__":
    main()
