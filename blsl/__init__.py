# Copyright (c) 2022 Dr. K. D. Murray/Gekkonid Consulting <spam@gekkonid.com> 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from sys import argv, exit, stderr

__version__ = "0.2.7"

cmds = {}


from .telogrep import telogrep_main
cmds["telogrep"] = telogrep_main

from .n50 import n50_main
cmds["n50"] = n50_main

from .falen import falen_main
cmds["falen"] = falen_main

from .mask2bed import mask2bed_main
cmds["mask2bed"] = mask2bed_main

from .pansn_rename import main as pansn_rename_main
cmds["pansn-rename"] = pansn_rename_main

from .genigvjs import genigvjs_main
cmds["genigvjs"] = genigvjs_main

from .ildemux import main as ildemux_main
cmds["ildemux"] = ildemux_main

from .ilsample import main as ilsample_main
cmds["ilsample"] = ilsample_main

from .regionbed import main as regionbed_main
cmds["regionbed"] = regionbed_main

from .uniref_accessionmap import main as uniref_main
cmds["uniref-acc2taxid"] = uniref_main

from .nstitch import nstitch_main
cmds["nstitch"] = nstitch_main

from .gg2k import gg2k_main
cmds["gg2k"] = gg2k_main

from .equalbestblast import equalbestblast_main
cmds["equalbestblast"] = equalbestblast_main

from .tabcat import tabcat_main
cmds["tabcat"] = tabcat_main

from .esearchandfetch import main as esf_main
cmds["esearchandfetch"] = esf_main

from .deepclust2fa import deepclust2fa_main
cmds["deepclust2fa"] = deepclust2fa_main

from .farename import farename_main
cmds["farename"] = farename_main

from .gffcat import gffcat_main
cmds["gffcat"] = gffcat_main

from .gffparse import gffparse_main
cmds["gffparse"] = gffparse_main

from .gffcsqify import main as gffcsqify_main
cmds["gffcsqify"] = gffcsqify_main

from .gfftagsane import main as gfftagsane_main
cmds["gfftagsane"] = gfftagsane_main

from .liftoff_gff3 import liftoff_gff3_main
cmds["liftoff-gff3"] = liftoff_gff3_main

from .ebiosra2rl2s import main as rl2s_main
cmds["ebiosra2rl2s"] = rl2s_main

from .galhist import main as galhist_main
cmds["galhist"] = galhist_main

from .pairslash import main as pairslash_main
cmds["pairslash"] = pairslash_main

from .pairs import main as pairs_main
cmds["pairs"] = pairs_main


try:
    from .vcfstats import main as vcfstats_main
    cmds["vcfstats"] = vcfstats_main

    from .vcfparallel import main as vcfparallel_main
    cmds["vcfparallel"] = vcfparallel_main
except ImportError as exc:
    if len(argv) < 2:
        print(str(exc), "-- disabling vcfstats command", file=stderr)

from .shannon import main as shannon_main
cmds["shannon-entropy"] = shannon_main

from .fastasanitiser import main as fastasanitiser_main
cmds["fastasanitiser"] = fastasanitiser_main

from .tidyqc import main as tidyqc_main
cmds["tidyqc"] = tidyqc_main

from .jsonl2csv import main as jsonl2csv_main
cmds["jsonl2csv"] = jsonl2csv_main


def mainhelp(argv=None):
    """Print this help message"""
    print("USAGE: blsl <subtool> [options...]\n\n")
    print("Where <subtool> is one of:\n")
    for tool, func in cmds.items():
        print("  {:<19}".format(tool + ":"), " ", func.__doc__.split("\n")[0])
    print("\n\nUse blsl subtool --help to get help about a specific tool")

cmds["help"] = mainhelp

def main():
    if len(argv) < 2:
        mainhelp()
        exit(0)
    if argv[1] not in cmds:
        print("ERROR:", argv[1], "is not a known subtool. See help below")
        mainhelp()
        exit(1)
    cmds[argv[1]](argv[2:])

