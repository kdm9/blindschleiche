# Copyright (c) 2022 Dr. K. D. Murray/Gekkonid Consulting <spam@gekkonid.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import argparse
import json

template = """
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


def genigvjs_main(argv=None):
    """Generate a simple IGV.js visualisation of some bioinf files."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--title", "-t", default="IGV.js",
            help="Webpage title")
    ap.add_argument("--reference", "-r",
            help="reference file. must be *.fa or *.fasta, with associated *.fai, and can also have *.gff as annotation.")
    ap.add_argument("--outdir", "-o",
            help="Output directory. Data files will be softlinked there, and 'index.html' will be generated.")
    ap.add_argument("tracks", help="Any files to be used as tracks. Can be gff/bed/vcf/bcf/bam/cram. Must be indexed", nargs="+")
    args = ap.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    data = {
        "reference": {},
        "tracks": [],
    }
    ref = Path(args.reference)
    refbase = ref.stem
    data["reference"] = {
        "id": refbase,
        "name": refbase,
        "fastaURL": f"./{ref.name}"
    }
    (outdir / ref.name).symlink_to(ref)
    (outdir / (ref.name + ".fai")).symlink_to(Path(str(ref) + ".fai"))

    
    for track in args.tracks:
        track = Path(track)
        base = track.stem
        (outdir / track.name).symlink_to(ref)
        trackdat = {
            "name": base,
            "url": f"./{track.name}",
            "format": track.suffix.lstrip("."),
        }
        data["tracks"].append(trackdat)

    with open(outdir / "index.html", "w") as fh:
        html = template \
                .replace("__title__", args.title) \
                .replace("__data__", json.dumps(data))
        fh.write(html)

if __name__ == "__main__":
    genigvjs_main()
