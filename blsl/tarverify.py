# Copyright (c) 2022 Dr. K. D. Murray/Gekkonid Consulting <spam@gekkonid.com>
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""Piece-wise tar archive verifier and builder.

Scans tar chunks in a directory, verifies each member against disk,
moves verified regular files and symlinks to trash, then archives
any remaining files as a new chunk.
"""

import argparse
import os
import re
import stat
import sys
import tarfile
import tempfile
import time


CHUNK_RE = re.compile(r'^chunk_(\d+)\.tar$')
CHUNK_FMT = 'chunk_{:04d}.tar'
COPY_BUFSIZE = 64 * 1024
PAX_FORMAT = tarfile.PAX_FORMAT


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description='Verify tar chunks against disk and archive remaining files.'
    )
    parser.add_argument(
        '-f', '--file', dest='tar_dir', required=True,
        help='Directory of tar chunk files'
    )
    parser.add_argument(
        '-t', '--trash', required=True,
        help='Trash directory for verified files'
    )
    parser.add_argument(
        '--timeout', type=float, metavar='SECS',
        help='Exit with code 2 after this many seconds'
    )
    parser.add_argument(
        'input', nargs='+',
        help='Input directories to archive'
    )
    return parser.parse_args(argv)


def resolve_disk_path(member_name, input_dirs):
    abs_path = '/' + member_name
    for d in input_dirs:
        nd = d.rstrip('/')
        if abs_path == nd or abs_path.startswith(nd + '/'):
            return abs_path
    return None


def stream_compare(fh_a, fh_b, bufsize=COPY_BUFSIZE):
    while True:
        a = fh_a.read(bufsize)
        b = fh_b.read(bufsize)
        if a != b:
            return False
        if not a:
            return True


def verify_member(tf, member, disk_path):
    try:
        lst = os.lstat(disk_path)
    except OSError:
        return False

    if member.isfile():
        if not stat.S_ISREG(lst.st_mode):
            return False
        if lst.st_size != member.size:
            return False
        if int(lst.st_mtime) != int(member.mtime):
            return False
        if member.size == 0:
            return True
        fh_tar = tf.extractfile(member)
        if fh_tar is None:
            return False
        try:
            with open(disk_path, 'rb') as fh_disk:
                return stream_compare(fh_tar, fh_disk)
        except OSError:
            return False

    if member.issym():
        if not stat.S_ISLNK(lst.st_mode):
            return False
        return os.readlink(disk_path) == member.linkname

    if member.isdir():
        return stat.S_ISDIR(lst.st_mode)

    return False


def move_to_trash(src, trash_dir):
    rel = src.lstrip('/')
    dst = os.path.join(trash_dir, rel)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    base = dst
    counter = 0
    while os.path.lexists(dst):
        counter += 1
        dst = f'{base}.{counter}'
    os.rename(src, dst)


def timeout_exceeded(timeout, start_time):
    if timeout is None:
        return False
    return time.monotonic() - start_time >= timeout


def next_chunk_number(tar_dir):
    max_num = -1
    if os.path.isdir(tar_dir):
        for entry in os.listdir(tar_dir):
            m = CHUNK_RE.match(entry)
            if m:
                num = int(m.group(1))
                if num > max_num:
                    max_num = num
    return max_num + 1


def get_chunks(tar_dir):
    n = -1
    if not os.path.isdir(tar_dir):
        return
    chunks = []
    for entry in os.listdir(tar_dir):
        if not entry.endswith(".tar"):
            continue
        m = CHUNK_RE.match(entry)
        if m:
            chunks.append((int(m.group(1)), entry))
        else:
            chunks.append((n, entry))
            n -= 1
    chunks.sort()
    for _, name in chunks:
        yield os.path.join(tar_dir, name)


def prune_excluded(dirpath, dirnames, excluded_real):
    to_remove = []
    for d in dirnames:
        full = os.path.join(dirpath, d)
        try:
            real = os.path.realpath(full)
        except OSError:
            continue
        if real in excluded_real:
            to_remove.append(d)
    for d in to_remove:
        dirnames.remove(d)


def phase1_verify(tar_dir, trash_dir, input_dirs, timeout, start_time):
    for chunk_path in get_chunks(tar_dir):
        try:
            tf = tarfile.open(chunk_path, 'r')
        except (tarfile.TarError, OSError):
            continue
        try:
            member = tf.next()
            while member:
                if timeout_exceeded(timeout, start_time):
                    return False

                disk_path = resolve_disk_path(member.name, input_dirs)
                if disk_path is None:
                    member = tf.next()
                    continue
                if not os.path.lexists(disk_path):
                    member = tf.next()
                    continue
                if verify_member(tf, member, disk_path):
                    if member.isfile() or member.issym():
                        print("VERIFY", member.name)
                        move_to_trash(disk_path, trash_dir)
                else:
                    print("FAIL_VERIFY", member.name)
                member = tf.next()
        finally:
            tf.close()
    return True


def phase2_append(tar_dir, trash_dir, input_dirs, timeout, start_time):
    excluded_real = set()
    for d in (trash_dir, tar_dir):
        try:
            excluded_real.add(os.path.realpath(d))
        except OSError:
            pass

    chunk_num = next_chunk_number(tar_dir)
    chunk_path = os.path.join(tar_dir, CHUNK_FMT.format(chunk_num))
    os.makedirs(tar_dir, exist_ok=True)

    fd, tmp_path = tempfile.mkstemp(dir=tar_dir, prefix='.tmp_chunk_', suffix='.tar')
    os.close(fd)

    any_added = False
    tf = None
    try:
        tf = tarfile.open(tmp_path, 'w', format=PAX_FORMAT)
        for input_dir in input_dirs:
            if not os.path.isdir(input_dir):
                continue
            for dirpath, dirnames, filenames in os.walk(input_dir):
                prune_excluded(dirpath, dirnames, excluded_real)

                for d in list(dirnames):
                    full_d = os.path.join(dirpath, d)
                    if os.path.islink(full_d):
                        dirnames.remove(d)
                        arcname = full_d.lstrip('/')
                        try:
                            tf.add(full_d, arcname=arcname, recursive=False)
                            any_added = True
                            print("ADD", arcname)
                        except (OSError, tarfile.TarError):
                            pass

                for f in filenames:
                    if timeout_exceeded(timeout, start_time):
                        return (False, False)
                    full_f = os.path.join(dirpath, f)
                    arcname = full_f.lstrip('/')
                    try:
                        if os.path.islink(full_f) or os.path.isfile(full_f):
                            tf.add(full_f, arcname=arcname)
                            any_added = True
                            print("ADD", arcname)
                    except (OSError, tarfile.TarError):
                        pass
    finally:
        if tf is not None:
            tf.close()
        if any_added:
            os.rename(tmp_path, chunk_path)
        else:
            os.unlink(tmp_path)

    return (True, any_added)


def main(argv=None):
    """Verify tar chunks against disk and archive remaining files."""
    args = parse_args(argv)
    start_time = time.monotonic()

    tar_dir = os.path.abspath(args.tar_dir)
    trash_dir = os.path.abspath(args.trash)
    input_dirs = [os.path.abspath(d) for d in args.input]
    timeout = args.timeout

    for d in input_dirs:
        if not os.path.isdir(d):
            print(f'error: input directory does not exist: {d}', file=sys.stderr)
            sys.exit(1)

    input_real = set()
    for d in input_dirs:
        try:
            input_real.add(os.path.realpath(d))
        except OSError:
            pass

    for name, d in (('trash', trash_dir), ('tar', tar_dir)):
        try:
            real_d = os.path.realpath(d)
            for ir in input_real:
                if real_d == ir or real_d.startswith(ir + '/'):
                    print(f'warning: {name} directory {d} is inside input '
                          f'directory {ir}', file=sys.stderr)
                    break
        except OSError:
            pass

    os.makedirs(trash_dir, exist_ok=True)

    ok1 = phase1_verify(tar_dir, trash_dir, input_dirs, timeout, start_time)
    if not ok1:
        sys.exit(2)

    ok2, any_added = phase2_append(tar_dir, trash_dir, input_dirs, timeout, start_time)
    if not ok2:
        sys.exit(2)
    if any_added:
        sys.exit(3)

    sys.exit(0)
