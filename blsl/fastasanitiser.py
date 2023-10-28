#!/usr/bin/env python3
# Copyright 2020 K.D. Murray
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


import re
from Bio import SeqIO
from hashlib import md5
import json


class FastaSanitiser(object):
    def __init__(self, regex=r'(S\d+)', forgiving=True, usemd5=False, old_idmap=None):
        self.regex = regex
        if isinstance(regex, str):
            self.regex = re.compile(regex)
        self.id_map = {}
        self.forgiving = forgiving
        self.usemd5 = usemd5

        if isinstance(old_idmap, str):
            with open(old_idmap) as fh:
                self.id_map = json.load(fh)
        elif isinstance(old_idmap, dict):
            self.id_map = old_idmap


    def _get_id(self, seqid):
        if self.usemd5:
            return md5(seqid.encode("utf-8")).hexdigest()
        m = self.regex.search(seqid)
        if m is not None:
            return m.group(0)
        elif m is None and self.forgiving:
            return md5(seqid.encode("utf-8")).hexdigest()
        else:
            raise ValueError(f"Couldn't match regex to sanitise name '{seqid}'")

    def sanitise(self, seqfile, format='fasta'):
        seqs = []
        with open(seqfile) as fh:
            for seq in SeqIO.parse(seqfile, format):
                print(seq.name, seq.id, seq.description)
                oldname = seq.description  # stupidly, some "ids" have spaces in them in gisaid data so we use the whole line.
                newname = self._get_id(oldname)
                if newname in self.id_map:
                    raise ValueError(f"Duplicate sequences with sanitised id {newname}")
                self.id_map[newname] = oldname
                seq.name = newname
                seq.id = newname
                seq.description = ""
                seqs.append(seq)
        return seqs

    def convert_forward(self, source, destination):
        conv_seqs = self.sanitise(source)
        with open(destination, "w") as seqfh, open(destination + ".idmap.json", "w") as mapfh:
            SeqIO.write(conv_seqs, seqfh, "fasta")
            json.dump(self.id_map, mapfh)

    def convert_back(self, source, destination):
        with open(source) as ifh, open(destination, "w") as ofh:
            for line in ifh:
                for id_new, id_old in self.id_map.items():
                    if id_new in line:
                        line = line.replace(id_new, id_old)
                        print(line, end="")
                ofh.write(line)


def main(argv=None):
    """Sanitise fasta IDs to something sane, then back again"""
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--backwards", "-b", action="store_true",
                    help="Convert from sanitised to unsantised format")
    ap.add_argument("--old-mapping", "-m", type=str,
                    help="Json file mapping clean IDs to old IDs.")
    ap.add_argument("--regex", "-r", type=str, default=r'(EPI_ISL_\d+)',
                    help="ID regex. Must capture the new ID as the first group.")
    ap.add_argument("--forgiving", "-f", action='store_true',
                    help="If --regex isn't found, use md5sum of ID instead.")
    ap.add_argument("--md5sum", action='store_true',
                    help="Instead of using --regex, always use md5sum of SeqIDs.")
    ap.add_argument("--output", "-o", type=str,
                    help="Output file")
    ap.add_argument("input", type=str,
                    help="Input file")
    args = ap.parse_args(argv)
    fs = FastaSanitiser(regex=args.regex, forgiving=args.forgiving, usemd5=args.md5sum, old_idmap=args.old_mapping)
    if not args.backwards:
        fs.convert_forward(args.input, args.output)
    else:
        fs.convert_back(args.input, args.output)


if __name__ == "__main__":
    main()
