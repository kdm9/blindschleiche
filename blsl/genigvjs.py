# Copyright (c) 2022 Dr. K. D. Murray/Gekkonid Consulting <spam@gekkonid.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import argparse
import json
from pathlib import Path
import subprocess
from base64 import b64encode 
import shutil


def get_data_uri(data):
    if isinstance(data, str) or isinstance(data, Path):
        data = Path(data)
        if data.stat().st_size > 2**22:
            raise ValueError(f"Too big to embed: {data}")
        with open(data, "rb") as fh:
            data = fh.read()
    elif len(data) > 2**22:
        raise ValueError(f"Too big to embed data")
    if data[0] == 0x1f and data[1] == 0x8b:
        mediatype = "data:application/gzip"
    else:
        mediatype = "data:application:octet-stream"

    enc_str = b64encode(data)

    data_uri = mediatype + ";base64," + str(enc_str)[2:-1]
    return data_uri


default_template = """
<!doctype html>
<html lang=en>
  <head>
    <meta charset=utf-8>
    <meta name=viewport content="width=device-width,initial-scale=1">
    <meta http-equiv=x-ua-compatible content="IE=edge,chrome=1">
    <title>__title__</title>
    <script src="https://cdn.jsdelivr.net/npm/igv@2.12.6/dist/igv.min.js"></script>
</head>
<body>
  <div id="igvherepls"></div>
  <script>
    var igvDiv = document.getElementById("igvherepls");
    var options = __data__;
        igv.createBrowser(igvDiv, options)
                .then(function (browser) {
                    console.log("Created IGV browser");
                })
  </script>
</body>
"""

def link(linkpath, target, copy=False):
    linkpath = Path(linkpath)
    target = Path(target).absolute()
    if copy:
        linkpath.parent.mkdir(exist_ok=True)
        shutil.copyfile(target, linkpath)
    else:
        try:
            if linkpath.samefile(target):
                return
        except:
            pass
        linkpath.symlink_to(target)

def genigvjs_main(argv=None):
    """Generate a simple IGV.js visualisation of some bioinf files."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", "-T", required=False,
            help="Alternative HTML template")
    ap.add_argument("--copy", "--cp", action="store_true",
            help="Copy, don't link, data to output dir")
    ap.add_argument("--embed", "-e", action="store_true",
            help="Encode data as base64 within html (almost always a bad idea on non-toy datasets).")
    ap.add_argument("--prefix", "-p", default="./",
            help="URL prefix of paths relative to generated index.html")
    ap.add_argument("--title", "-t", default="IGV.js",
            help="Webpage title")
    ap.add_argument("--reference", "-r", required=True,
            help="reference file. must be *.fa or *.fasta, with associated *.fai, and can also have *.gff as annotation.")
    ap.add_argument("--locus", "-l", 
            help="Initial region of interest.")
    ap.add_argument("--outdir", "-o", required=True,
            help="Output directory. Data files will be softlinked there, and 'index.html' will be generated.")
    ap.add_argument("tracks", help="Any files to be used as tracks. Can be gff/bed/vcf/bcf/bam/cram. Must be indexed", nargs="+")
    args = ap.parse_args(argv)

    if args.template is not None:
        with open(args.template) as fh:
            template = fh.read()
    else:
        template = default_template

    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    data = {
        "reference": {},
        "tracks": [],
    }
    if args.locus:
        data["locus"] = args.locus
    ref = Path(args.reference)
    refbase = ref.stem
    reffai = Path(str(ref) + ".fai")
    data["reference"] = {
        "id": refbase,
        "name": refbase,
    }
    if args.embed:
        data["reference"].update({
            "fastaURL": get_data_uri(ref),
            "indexURL": get_data_uri(reffai),
        })
    else:
        data["reference"].update({
            "fastaURL": f"{args.prefix}{ref.name}",
        })
        link(outdir / ref.name, ref, copy=args.copy)
        link(outdir / (ref.name + ".fai"), reffai, copy=args.copy)
    cbpaired = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a', '#b15928']
    for i, track in enumerate(args.tracks):
        track = Path(track)
        base = track.stem
        dots = track.name.split('.')
        format = dots[-1]
        if format == "gz":
            format = dots[-2]
        index = None
        if format == "bam":
            index = Path(str(track) + ".bai")
        else:
            index = Path(str(track) + ".tbi")
        if not index.exists():
            index = None
        trackdat = {
            "name": base,
            "format": format,
            "autoHeight": True,
            "minHeight": 50,
            "maxHeight": 500,
            "color": cbpaired[i%len(cbpaired)],
        }
        if args.embed:
            trackdat["url"] = get_data_uri(track)
            if index is not None:
                trackdat["indexURL"] = get_data_uri(index)
        else:
            trackdat["url"] =  f"{args.prefix}{track.name}"
            link(outdir / track.name, track, copy=args.copy)
            if index is not None:
                trackdat["indexURL"] =  f"{args.prefix}{index.name}"
                link(outdir / index.name, index, copy=args.copy)
        data["tracks"].append(trackdat)
    with open(outdir / "index.html", "w") as fh:
        html = template \
                .replace("__title__", args.title) \
                .replace("__data__", json.dumps(data, indent=4))
        fh.write(html)

if __name__ == "__main__":
    genigvjs_main()
