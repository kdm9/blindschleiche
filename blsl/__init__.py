# Copyright (c) 2022 Dr. K. D. Murray/Gekkonid Consulting <spam@gekkonid.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from sys import argv, exit

__version__ = "0.1.7"

cmds = {}

from .telogrep import telogrep_main
cmds["telogrep"] = telogrep_main

from .n50 import n50_main
cmds["n50"] = n50_main

from .falen import falen_main
cmds["falen"] = falen_main

from .mask2bed import mask2bed_main
cmds["mask2bed"] = mask2bed_main

from .genigvjs import genigvjs_main
cmds["genigvjs"] = genigvjs_main

from .liftoff_gff3 import liftoff_gff3_main
cmds["liftoff-gff3"] = liftoff_gff3_main

from .pansn_rename import main as pansn_rename_main
cmds["pansn-rename"] = pansn_rename_main

from .ildemux import main as ildemux_main
cmds["ildemux"] = ildemux_main


def mainhelp():
    print("USAGE: blsl <subtool> [options...]\n\n")
    print("Where <subtool> is one of:\n")
    for tool, func in cmds.items():
        print("  {:<15}".format(tool + ":"), " ", func.__doc__.split("\n")[0])
    print("\n\nUse blsl subtool --help to get help about a specific tool")


def main():
    if len(argv) < 2:
        mainhelp()
        exit(0)
    if argv[1] not in cmds:
        print("ERROR:", argv[1], "is not a known subtool. See help below")
        mainhelp()
        exit(1)
    cmds[argv[1]](argv[2:])

