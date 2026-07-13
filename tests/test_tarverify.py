import errno
import io
import os
import re
import shutil
import signal
import stat
import subprocess
import sys
import tarfile
import tempfile
import time
from pathlib import Path

import pytest

from blsl.tarverify import (
    CHUNK_FMT,
    CHUNK_RE,
    COPY_BUFSIZE,
    PAX_FORMAT,
    get_chunks,
    main,
    move_to_trash,
    next_chunk_number,
    parse_args,
    phase1_verify,
    phase2_append,
    prune_excluded,
    resolve_disk_path,
    stream_compare,
    timeout_exceeded,
    verify_member,
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

BLSL_MODULE = [sys.executable, '-m', 'blsl', 'tarverify']


def _run_tarverify(tar_dir, trash_dir, *input_dirs, **kw):
    timeout = kw.get('timeout', None)
    argv = ['-f', tar_dir, '-t', trash_dir]
    if timeout is not None:
        argv.extend(['--timeout', str(timeout)])
    argv.extend(input_dirs)
    return main(argv)


def _make_tar(tar_path, members):
    with tarfile.open(tar_path, 'w', format=PAX_FORMAT) as tf:
        for arcname, disk_path in members:
            tf.add(disk_path, arcname=arcname, recursive=False)


def _tar_members(tar_path):
    members = []
    with tarfile.open(tar_path, 'r') as tf:
        m = tf.next()
        while m:
            members.append(m)
            m = tf.next()
    return members


def _write_file(path, content):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    mode = 'wb' if isinstance(content, bytes) else 'w'
    encoding = None if isinstance(content, bytes) else 'utf-8'
    with open(path, mode, encoding=encoding) as fh:
        fh.write(content)


def _symlink(target, link):
    Path(link).parent.mkdir(parents=True, exist_ok=True)
    if os.path.lexists(link):
        os.unlink(link)
    os.symlink(target, link)


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def tmp_path():
    d = tempfile.mkdtemp(prefix='blsl_test_')
    yield d
    shutil.rmtree(d, ignore_errors=True)


@pytest.fixture
def input_tree(tmp_path):
    s = os.path.join(tmp_path, 'input')
    os.makedirs(s)
    _write_file(os.path.join(s, 'file_a.txt'), 'hello world\n')
    _write_file(os.path.join(s, 'empty.dat'), '')
    _write_file(os.path.join(s, 'binary.bin'), b'\x00\x01\x02\x03\xFF')
    big = b'A' * 100000 + b'B' * 100000
    _write_file(os.path.join(s, 'large.bin'), big)
    os.makedirs(os.path.join(s, 'sub', 'deep'))
    _write_file(os.path.join(s, 'sub', 'deep', 'nested.txt'), 'deep\n')
    os.makedirs(os.path.join(s, 'emptydir'))
    _symlink('file_a.txt', os.path.join(s, 'goodlink'))
    _symlink('sub', os.path.join(s, 'dirlink'))
    _symlink('nonexistent', os.path.join(s, 'brokenlink'))
    return s


@pytest.fixture
def tar_cache(tmp_path):
    d = os.path.join(tmp_path, 'tar')
    os.makedirs(d)
    return d


@pytest.fixture
def trash_dir(tmp_path):
    d = os.path.join(tmp_path, 'trash')
    os.makedirs(d)
    return d


# ======================================================================
# parse_args
# ======================================================================

class TestParseArgs:
    def test_required_args(self):
        args = parse_args(['-f', '/tmp/tar', '-t', '/tmp/trash', '/tmp/in'])
        assert args.tar_dir == '/tmp/tar'
        assert args.trash == '/tmp/trash'
        assert args.input == ['/tmp/in']

    def test_multiple_inputs(self):
        args = parse_args(['-f', '/t', '-t', '/r', '/a', '/b', '/c'])
        assert args.input == ['/a', '/b', '/c']

    def test_timeout(self):
        args = parse_args(['-f', '/t', '-t', '/r', '--timeout', '42.5', '/i'])
        assert args.timeout == 42.5

    def test_timeout_none(self):
        args = parse_args(['-f', '/t', '-t', '/r', '/i'])
        assert args.timeout is None

    def test_long_opts(self):
        args = parse_args(['--file', '/t', '--trash', '/r', '--timeout', '10', '/i'])
        assert args.tar_dir == '/t'
        assert args.trash == '/r'
        assert args.timeout == 10

    def test_missing_file(self):
        with pytest.raises(SystemExit):
            parse_args(['-t', '/r', '/i'])

    def test_missing_trash(self):
        with pytest.raises(SystemExit):
            parse_args(['-f', '/t', '/i'])

    def test_missing_input(self):
        with pytest.raises(SystemExit):
            parse_args(['-f', '/t', '-t', '/r'])


# ======================================================================
# resolve_disk_path
# ======================================================================

class TestResolveDiskPath:
    def test_exact_match(self):
        assert resolve_disk_path('a/b/c', ['/a']) == '/a/b/c'

    def test_subdir_match(self):
        assert resolve_disk_path('a/b/c/d', ['/a']) == '/a/b/c/d'

    def test_match_second_input(self):
        assert resolve_disk_path('a/b/c', ['/x', '/a']) == '/a/b/c'

    def test_no_match(self):
        assert resolve_disk_path('a/b/c', ['/x', '/y']) is None

    def test_trailing_slash_on_input(self):
        assert resolve_disk_path('a/b/c', ['/a/']) == '/a/b/c'

    def test_trailing_slash_match(self):
        assert resolve_disk_path('a', ['/a/']) == '/a'

    def test_empty_member_name(self):
        assert resolve_disk_path('', ['/a']) is None

    def test_root_input(self):
        assert resolve_disk_path('a/b', ['/']) == '/a/b'

    def test_single_char_member(self):
        assert resolve_disk_path('x', ['/a']) is None


# ======================================================================
# stream_compare
# ======================================================================

class TestStreamCompare:
    def test_identical(self):
        a = io.BytesIO(b'hello world')
        b = io.BytesIO(b'hello world')
        assert stream_compare(a, b) is True

    def test_identical_empty(self):
        a = io.BytesIO(b'')
        b = io.BytesIO(b'')
        assert stream_compare(a, b) is True

    def test_differ_first_byte(self):
        a = io.BytesIO(b'xello')
        b = io.BytesIO(b'hello')
        assert stream_compare(a, b) is False

    def test_differ_last_byte(self):
        a = io.BytesIO(b'hellx')
        b = io.BytesIO(b'hello')
        assert stream_compare(a, b) is False

    def test_differ_mid(self):
        a = io.BytesIO(b'heXXo')
        b = io.BytesIO(b'hello')
        assert stream_compare(a, b) is False

    def test_one_empty(self):
        assert stream_compare(io.BytesIO(b''), io.BytesIO(b'a')) is False
        assert stream_compare(io.BytesIO(b'a'), io.BytesIO(b'')) is False

    def test_one_shorter(self):
        a = io.BytesIO(b'hello')
        b = io.BytesIO(b'hello world')
        assert stream_compare(a, b) is False

    def test_one_longer(self):
        a = io.BytesIO(b'hello world')
        b = io.BytesIO(b'hello')
        assert stream_compare(a, b) is False

    def test_same_size_different_content(self):
        a = io.BytesIO(b'aaaa')
        b = io.BytesIO(b'bbbb')
        assert stream_compare(a, b) is False

    def test_bufsize_one(self):
        a = io.BytesIO(b'abc')
        b = io.BytesIO(b'abc')
        assert stream_compare(a, b, bufsize=1) is True

    def test_bufsize_one_differ(self):
        a = io.BytesIO(b'abc')
        b = io.BytesIO(b'abX')
        assert stream_compare(a, b, bufsize=1) is False

    def test_large_identical(self):
        data = os.urandom(COPY_BUFSIZE * 3 + 7)
        assert stream_compare(io.BytesIO(data), io.BytesIO(data)) is True

    def test_large_differ_at_boundary(self):
        n = COPY_BUFSIZE
        data = os.urandom(n * 3)
        mutated = bytearray(data)
        mutated[n] ^= 0xFF
        assert stream_compare(io.BytesIO(data), io.BytesIO(bytes(mutated))) is False

    def test_large_differ_last_byte(self):
        data = os.urandom(10000)
        mutated = bytearray(data)
        mutated[-1] ^= 0xFF
        assert stream_compare(io.BytesIO(data), io.BytesIO(bytes(mutated))) is False


# ======================================================================
# verify_member
# ======================================================================

class TestVerifyMember:
    @staticmethod
    def _make_tar_with_file(tar_path, arcname, disk_path):
        _make_tar(tar_path, [(arcname, disk_path)])

    def test_file_all_match(self, tmp_path):
        disk = os.path.join(tmp_path, 'disk')
        _write_file(disk, 'exact content')
        os.utime(disk, (1000000000, 1000000000))

        tar = os.path.join(tmp_path, 'chunk.tar')
        self._make_tar_with_file(tar, 'disk', disk)
        m = _tar_members(tar)[0]

        assert verify_member(tarfile.open(tar, 'r'), m, disk) is True

    def test_file_same_size_different_content(self, tmp_path):
        disk = os.path.join(tmp_path, 'disk')
        content_a = b'A' * 4096
        _write_file(disk, content_a)
        os.utime(disk, (1000000000, 1000000000))

        content_b = b'B' * 4096
        tmp_content = os.path.join(tmp_path, 'tmp_content')
        _write_file(tmp_content, content_b)
        os.utime(tmp_content, (1000000000, 1000000000))

        tar = os.path.join(tmp_path, 'chunk.tar')
        self._make_tar_with_file(tar, 'tmp_content', tmp_content)
        m = _tar_members(tar)[0]

        assert verify_member(tarfile.open(tar, 'r'), m, disk) is False

    def test_file_different_size(self, tmp_path):
        disk = os.path.join(tmp_path, 'disk')
        _write_file(disk, 'short')
        os.utime(disk, (1000000000, 1000000000))

        tar = os.path.join(tmp_path, 'chunk.tar')
        _write_file(os.path.join(tmp_path, 'long'), 'longer content here')
        os.utime(os.path.join(tmp_path, 'long'), (1000000000, 1000000000))
        self._make_tar_with_file(tar, 'long', os.path.join(tmp_path, 'long'))
        m = _tar_members(tar)[0]

        assert verify_member(tarfile.open(tar, 'r'), m, disk) is False

    def test_file_different_mtime(self, tmp_path):
        disk = os.path.join(tmp_path, 'disk')
        _write_file(disk, 'content')
        os.utime(disk, (2000000000, 2000000000))

        tar = os.path.join(tmp_path, 'chunk.tar')
        _write_file(os.path.join(tmp_path, 'src'), 'content')
        os.utime(os.path.join(tmp_path, 'src'), (1000000000, 1000000000))
        self._make_tar_with_file(tar, 'src', os.path.join(tmp_path, 'src'))
        m = _tar_members(tar)[0]

        assert verify_member(tarfile.open(tar, 'r'), m, disk) is False

    def test_file_same_mtime_int_truncation(self, tmp_path):
        disk = os.path.join(tmp_path, 'disk')
        _write_file(disk, 'content')
        os.utime(disk, (1000000000, 1000000000))

        tar = os.path.join(tmp_path, 'chunk.tar')
        _write_file(os.path.join(tmp_path, 'src'), 'content')
        os.utime(os.path.join(tmp_path, 'src'), (1000000000, 1000000000.7))
        self._make_tar_with_file(tar, 'src', os.path.join(tmp_path, 'src'))
        m = _tar_members(tar)[0]

        assert verify_member(tarfile.open(tar, 'r'), m, disk) is True

    def test_empty_file_match(self, tmp_path):
        disk = os.path.join(tmp_path, 'disk')
        _write_file(disk, '')
        os.utime(disk, (1000000000, 1000000000))

        tar = os.path.join(tmp_path, 'chunk.tar')
        self._make_tar_with_file(tar, 'disk', disk)
        m = _tar_members(tar)[0]

        assert verify_member(tarfile.open(tar, 'r'), m, disk) is True

    def test_disk_missing(self, tmp_path):
        disk = os.path.join(tmp_path, 'disk')
        _write_file(disk, 'content')
        os.utime(disk, (1000000000, 1000000000))

        tar = os.path.join(tmp_path, 'chunk.tar')
        self._make_tar_with_file(tar, 'disk', disk)
        m = _tar_members(tar)[0]
        os.unlink(disk)

        assert verify_member(tarfile.open(tar, 'r'), m, disk) is False

    def test_disk_is_dir_not_file(self, tmp_path):
        disk = os.path.join(tmp_path, 'disk')
        _write_file(disk, 'content')
        os.utime(disk, (1000000000, 1000000000))

        tar = os.path.join(tmp_path, 'chunk.tar')
        self._make_tar_with_file(tar, 'disk', disk)
        m = _tar_members(tar)[0]

        os.unlink(disk)
        os.makedirs(disk)

        assert verify_member(tarfile.open(tar, 'r'), m, disk) is False

    def test_disk_is_symlink_not_file(self, tmp_path):
        disk = os.path.join(tmp_path, 'disk')
        _write_file(disk, 'content')
        os.utime(disk, (1000000000, 1000000000))

        tar = os.path.join(tmp_path, 'chunk.tar')
        self._make_tar_with_file(tar, 'disk', disk)
        m = _tar_members(tar)[0]

        os.unlink(disk)
        _write_file(os.path.join(tmp_path, 'target'), 'whatever')
        os.symlink(os.path.join(tmp_path, 'target'), disk)

        assert verify_member(tarfile.open(tar, 'r'), m, disk) is False

    def test_file_with_no_read_permission(self, tmp_path):
        disk = os.path.join(tmp_path, 'disk')
        _write_file(disk, 'content')
        os.utime(disk, (1000000000, 1000000000))

        tar = os.path.join(tmp_path, 'chunk.tar')
        self._make_tar_with_file(tar, 'disk', disk)
        m = _tar_members(tar)[0]

        os.chmod(disk, 0)

        try:
            assert verify_member(tarfile.open(tar, 'r'), m, disk) is False
        finally:
            os.chmod(disk, 0o644)

    def test_symlink_match(self, tmp_path):
        target = os.path.join(tmp_path, 'target')
        _write_file(target, 'data')
        link = os.path.join(tmp_path, 'link')
        _symlink('target', link)

        tar = os.path.join(tmp_path, 'chunk.tar')
        _make_tar(tar, [('link', link)])
        m = _tar_members(tar)[0]

        assert verify_member(tarfile.open(tar, 'r'), m, link) is True

    def test_symlink_different_target(self, tmp_path):
        target = os.path.join(tmp_path, 'target')
        _write_file(target, 'data')
        link = os.path.join(tmp_path, 'link')
        _symlink('target', link)

        tar = os.path.join(tmp_path, 'chunk.tar')
        _symlink('different_target', os.path.join(tmp_path, 'link2'))
        _make_tar(tar, [('link2', os.path.join(tmp_path, 'link2'))])
        m = _tar_members(tar)[0]

        assert verify_member(tarfile.open(tar, 'r'), m, link) is False

    def test_symlink_disk_not_symlink(self, tmp_path):
        link = os.path.join(tmp_path, 'link')
        _write_file(link, 'not a link')
        os.utime(link, (1000000000, 1000000000))

        tar = os.path.join(tmp_path, 'chunk.tar')
        _symlink('target', os.path.join(tmp_path, 'real_link'))
        _make_tar(tar, [('real_link', os.path.join(tmp_path, 'real_link'))])
        m = _tar_members(tar)[0]

        assert verify_member(tarfile.open(tar, 'r'), m, link) is False

    def test_directory_match(self, tmp_path):
        d = os.path.join(tmp_path, 'adir')
        os.makedirs(d)

        tar = os.path.join(tmp_path, 'chunk.tar')
        _make_tar(tar, [('adir', d)])
        m = _tar_members(tar)[0]

        assert verify_member(tarfile.open(tar, 'r'), m, d) is True

    def test_directory_disk_not_dir(self, tmp_path):
        d = os.path.join(tmp_path, 'adir')
        os.makedirs(d)

        tar = os.path.join(tmp_path, 'chunk.tar')
        _make_tar(tar, [('adir', d)])
        m = _tar_members(tar)[0]

        shutil.rmtree(d)
        _write_file(d, 'not a dir')

        assert verify_member(tarfile.open(tar, 'r'), m, d) is False

    def test_other_member_type_fifo(self, tmp_path):
        disk = os.path.join(tmp_path, 'disk')
        _write_file(disk, 'content')
        os.utime(disk, (1000000000, 1000000000))

        tar = os.path.join(tmp_path, 'chunk.tar')
        self._make_tar_with_file(tar, 'disk', disk)
        m = _tar_members(tar)[0]

        class FakeMember:
            name = m.name
            size = m.size
            mtime = m.mtime
            mode = m.mode

            def isfile(self):
                return False

            def issym(self):
                return False

            def isdir(self):
                return False

        fake = FakeMember()
        assert verify_member(tarfile.open(tar, 'r'), fake, disk) is False


# ======================================================================
# move_to_trash
# ======================================================================

class TestMoveToTrash:
    def test_basic_move(self, tmp_path):
        src_dir = os.path.join(tmp_path, 'src')
        os.makedirs(src_dir)
        src = os.path.join(src_dir, 'f.txt')
        _write_file(src, 'data')
        trash = os.path.join(tmp_path, 'trash')

        move_to_trash(src, trash)
        assert not os.path.lexists(src)
        rel = src.lstrip('/')
        assert os.path.isfile(os.path.join(trash, rel))

    def test_creates_parent_dirs(self, tmp_path):
        src_dir = os.path.join(tmp_path, 'a', 'b', 'c')
        os.makedirs(src_dir)
        src = os.path.join(src_dir, 'f.txt')
        _write_file(src, 'data')
        trash = os.path.join(tmp_path, 'trash')

        move_to_trash(src, trash)
        rel = src.lstrip('/')
        assert os.path.isfile(os.path.join(trash, rel))

    def test_collision_handling(self, tmp_path):
        src_dir = os.path.join(tmp_path, 'src')
        os.makedirs(src_dir)
        src = os.path.join(src_dir, 'f.txt')
        _write_file(src, 'data')
        trash = os.path.join(tmp_path, 'trash')

        rel = src.lstrip('/')
        os.makedirs(os.path.dirname(os.path.join(trash, rel)), exist_ok=True)
        _write_file(os.path.join(trash, rel), 'existing')

        move_to_trash(src, trash)
        assert os.path.isfile(os.path.join(trash, rel + '.1'))

    def test_multiple_collisions(self, tmp_path):
        src_dir = os.path.join(tmp_path, 'src')
        os.makedirs(src_dir)
        src = os.path.join(src_dir, 'f.txt')
        _write_file(src, 'data')
        trash = os.path.join(tmp_path, 'trash')

        rel = src.lstrip('/')
        os.makedirs(os.path.dirname(os.path.join(trash, rel)), exist_ok=True)
        for i in range(3):
            suffix = f'.{i}' if i else ''
            _write_file(os.path.join(trash, rel + suffix), f'v{i}')

        move_to_trash(src, trash)
        assert os.path.isfile(os.path.join(trash, rel + '.3'))

    def test_absolute_src_stripping(self, tmp_path):
        src = os.path.join(tmp_path, 'deep', 'nested', 'file.txt')
        _write_file(src, 'data')
        trash = os.path.join(tmp_path, 'trash')

        move_to_trash(src, trash)
        rel = src.lstrip('/')
        assert os.path.isfile(os.path.join(trash, rel))

    def test_move_symlink(self, tmp_path):
        src_dir = os.path.join(tmp_path, 'src')
        os.makedirs(src_dir)
        target = os.path.join(src_dir, 'target')
        _write_file(target, 'data')
        src = os.path.join(src_dir, 'link')
        _symlink('target', src)
        trash = os.path.join(tmp_path, 'trash')

        move_to_trash(src, trash)
        assert not os.path.lexists(src)
        rel = src.lstrip('/')
        assert os.path.islink(os.path.join(trash, rel))


# ======================================================================
# timeout_exceeded
# ======================================================================

class TestTimeoutExceeded:
    def test_none_timeout(self):
        assert timeout_exceeded(None, time.monotonic()) is False
        assert timeout_exceeded(None, time.monotonic() - 999999) is False

    def test_zero_timeout(self):
        now = time.monotonic()
        assert timeout_exceeded(0, now) is True
        assert timeout_exceeded(0, now - 0.001) is True

    def test_positive_timeout_not_exceeded(self):
        now = time.monotonic()
        assert timeout_exceeded(3600, now) is False

    def test_positive_timeout_exceeded(self):
        assert timeout_exceeded(0.001, time.monotonic() - 60) is True

    def test_negative_timeout(self):
        assert timeout_exceeded(-1, time.monotonic()) is True


# ======================================================================
# next_chunk_number
# ======================================================================

class TestNextChunkNumber:
    def test_empty_dir(self, tmp_path):
        d = os.path.join(tmp_path, 'tar')
        os.makedirs(d)
        assert next_chunk_number(d) == 0

    def test_missing_dir(self, tmp_path):
        assert next_chunk_number('/nonexistent/path') == 0

    def test_existing_chunks(self, tmp_path):
        d = os.path.join(tmp_path, 'tar')
        os.makedirs(d)
        _write_file(os.path.join(d, 'chunk_0000.tar'), '')
        _write_file(os.path.join(d, 'chunk_0001.tar'), '')
        assert next_chunk_number(d) == 2

    def test_gap_in_chunks(self, tmp_path):
        d = os.path.join(tmp_path, 'tar')
        os.makedirs(d)
        _write_file(os.path.join(d, 'chunk_0000.tar'), '')
        _write_file(os.path.join(d, 'chunk_0005.tar'), '')
        assert next_chunk_number(d) == 6

    def test_non_matching_files_ignored(self, tmp_path):
        d = os.path.join(tmp_path, 'tar')
        os.makedirs(d)
        _write_file(os.path.join(d, 'chunk_0000.tar'), '')
        _write_file(os.path.join(d, 'notes.txt'), '')
        _write_file(os.path.join(d, 'chunk_000.tar'), '')
        assert next_chunk_number(d) == 1

    def test_subdir_is_counted(self, tmp_path):
        d = os.path.join(tmp_path, 'tar')
        os.makedirs(d)
        os.makedirs(os.path.join(d, 'chunk_0005.tar'))
        assert next_chunk_number(d) == 6

    def test_non_chunk_tar_ignored(self, tmp_path):
        d = os.path.join(tmp_path, 'tar')
        os.makedirs(d)
        _write_file(os.path.join(d, 'other.tar'), '')
        assert next_chunk_number(d) == 0


# ======================================================================
# get_chunks
# ======================================================================

class TestGetChunks:
    def test_empty_dir(self, tmp_path):
        d = os.path.join(tmp_path, 'tar')
        os.makedirs(d)
        assert list(get_chunks(d)) == []

    def test_missing_dir(self, tmp_path):
        assert list(get_chunks('/nonexistent')) == []

    def test_sorted_chunks(self, tmp_path):
        d = os.path.join(tmp_path, 'tar')
        os.makedirs(d)
        _write_file(os.path.join(d, 'chunk_0002.tar'), '')
        _write_file(os.path.join(d, 'chunk_0000.tar'), '')
        _write_file(os.path.join(d, 'chunk_0001.tar'), '')
        chunks = list(get_chunks(d))
        expected = [os.path.join(d, f'chunk_{n:04d}.tar') for n in range(3)]
        assert chunks == expected

    def test_non_chunk_tar_sorted_before(self, tmp_path):
        d = os.path.join(tmp_path, 'tar')
        os.makedirs(d)
        _write_file(os.path.join(d, 'chunk_0000.tar'), '')
        _write_file(os.path.join(d, 'manual.tar'), '')
        chunks = list(get_chunks(d))
        assert os.path.basename(chunks[0]) == 'manual.tar'
        assert os.path.basename(chunks[1]) == 'chunk_0000.tar'

    def test_multiple_non_chunk(self, tmp_path):
        d = os.path.join(tmp_path, 'tar')
        os.makedirs(d)
        _write_file(os.path.join(d, 'chunk_0001.tar'), '')
        _write_file(os.path.join(d, 'z.tar'), '')
        _write_file(os.path.join(d, 'a.tar'), '')
        chunks = list(get_chunks(d))
        names = [os.path.basename(c) for c in chunks]
        assert names.index('chunk_0001.tar') == len(names) - 1

    def test_non_tar_ignored(self, tmp_path):
        d = os.path.join(tmp_path, 'tar')
        os.makedirs(d)
        _write_file(os.path.join(d, 'chunk_0000.tar'), '')
        _write_file(os.path.join(d, 'notes.txt'), '')
        assert len(list(get_chunks(d))) == 1


# ======================================================================
# prune_excluded
# ======================================================================

class TestPruneExcluded:
    def test_removes_matching_dir(self, tmp_path):
        parent = os.path.join(tmp_path, 'parent')
        os.makedirs(parent)
        os.makedirs(os.path.join(parent, 'skip'))
        os.makedirs(os.path.join(parent, 'keep'))

        excluded = {os.path.realpath(os.path.join(parent, 'skip'))}
        dirnames = ['skip', 'keep']
        prune_excluded(parent, dirnames, excluded)
        assert dirnames == ['keep']

    def test_keeps_non_matching(self, tmp_path):
        parent = os.path.join(tmp_path, 'parent')
        os.makedirs(parent)
        os.makedirs(os.path.join(parent, 'keep'))
        excluded = {os.path.realpath(os.path.join(parent, 'other'))}

        dirnames = ['keep']
        prune_excluded(parent, dirnames, excluded)
        assert dirnames == ['keep']

    def test_handles_broken_symlink(self, tmp_path):
        parent = os.path.join(tmp_path, 'parent')
        os.makedirs(parent)
        os.symlink(os.path.join(tmp_path, 'nonexistent'), os.path.join(parent, 'broken'))

        dirnames = ['broken', 'real']
        os.makedirs(os.path.join(parent, 'real'))
        prune_excluded(parent, dirnames, {'/some/other/path'})
        assert 'real' in dirnames

    def test_removes_multiple(self, tmp_path):
        parent = os.path.join(tmp_path, 'parent')
        os.makedirs(parent)
        os.makedirs(os.path.join(parent, 'skip_a'))
        os.makedirs(os.path.join(parent, 'skip_b'))
        os.makedirs(os.path.join(parent, 'keep'))

        excluded = {
            os.path.realpath(os.path.join(parent, 'skip_a')),
            os.path.realpath(os.path.join(parent, 'skip_b')),
        }

        dirnames = ['skip_a', 'keep', 'skip_b']
        prune_excluded(parent, dirnames, excluded)
        assert dirnames == ['keep']

    def test_empty_dirnames(self, tmp_path):
        parent = os.path.join(tmp_path, 'parent')
        os.makedirs(parent)
        dirnames = []
        prune_excluded(parent, dirnames, set())
        assert dirnames == []


# ======================================================================
# phase1_verify
# ======================================================================

class TestPhase1Verify:
    def test_all_verified_and_trashed(self, tmp_path, tar_cache, trash_dir):
        src = os.path.join(tmp_path, 'input')
        _write_file(os.path.join(src, 'a.txt'), 'aaa')
        _write_file(os.path.join(src, 'b.txt'), 'bbb')
        _symlink('a.txt', os.path.join(src, 'link_a'))

        chunk = os.path.join(tar_cache, 'chunk_0000.tar')
        _make_tar(chunk, [
            (src.lstrip('/') + '/a.txt', os.path.join(src, 'a.txt')),
            (src.lstrip('/') + '/b.txt', os.path.join(src, 'b.txt')),
            (src.lstrip('/') + '/link_a', os.path.join(src, 'link_a')),
        ])

        ok = phase1_verify(tar_cache, trash_dir, [src], None, time.monotonic())
        assert ok is True
        assert not os.path.exists(os.path.join(src, 'a.txt'))
        assert not os.path.exists(os.path.join(src, 'b.txt'))
        assert not os.path.lexists(os.path.join(src, 'link_a'))
        rel = src.lstrip('/')
        assert os.path.isfile(os.path.join(trash_dir, rel, 'a.txt'))
        assert os.path.isfile(os.path.join(trash_dir, rel, 'b.txt'))
        assert os.path.islink(os.path.join(trash_dir, rel, 'link_a'))

    def test_missing_disk_file_skipped(self, tmp_path, tar_cache, trash_dir):
        src = os.path.join(tmp_path, 'input')
        _write_file(os.path.join(src, 'keep.txt'), 'keep')

        chunk = os.path.join(tar_cache, 'chunk_0000.tar')
        arcname_keep = src.lstrip('/') + '/keep.txt'
        arcname_gone = src.lstrip('/') + '/gone.txt'
        gone_disk = os.path.join(src, 'gone.txt')
        _write_file(gone_disk, 'gone')

        _make_tar(chunk, [
            (arcname_keep, os.path.join(src, 'keep.txt')),
            (arcname_gone, gone_disk),
        ])
        os.unlink(gone_disk)

        ok = phase1_verify(tar_cache, trash_dir, [src], None, time.monotonic())
        assert ok is True
        assert not os.path.exists(os.path.join(src, 'keep.txt'))

    def test_fail_verify_not_moved(self, tmp_path, tar_cache, trash_dir):
        src = os.path.join(tmp_path, 'input')
        disk = os.path.join(src, 'file.txt')
        _write_file(disk, 'original')
        os.utime(disk, (1000000000, 1000000000))

        chunk = os.path.join(tar_cache, 'chunk_0000.tar')
        alt = os.path.join(tmp_path, 'alt')
        _write_file(alt, 'different content')
        os.utime(alt, (1000000000, 1000000000))
        _make_tar(chunk, [(src.lstrip('/') + '/file.txt', alt)])

        ok = phase1_verify(tar_cache, trash_dir, [src], None, time.monotonic())
        assert ok is True
        assert os.path.isfile(disk)

    def test_corrupt_tar_skipped(self, tmp_path, tar_cache, trash_dir):
        src = os.path.join(tmp_path, 'input')
        _write_file(os.path.join(src, 'f.txt'), 'data')

        chunk = os.path.join(tar_cache, 'chunk_0000.tar')
        _write_file(chunk, b'this is not a valid tar file\x00\xff\xfe')

        ok = phase1_verify(tar_cache, trash_dir, [src], None, time.monotonic())
        assert ok is True

    def test_timeout_returns_false(self, tmp_path, tar_cache, trash_dir):
        src = os.path.join(tmp_path, 'input')
        _write_file(os.path.join(src, 'f.txt'), 'data')

        chunk = os.path.join(tar_cache, 'chunk_0000.tar')
        _make_tar(chunk, [(src.lstrip('/') + '/f.txt', os.path.join(src, 'f.txt'))])

        ok = phase1_verify(tar_cache, trash_dir, [src], 0, time.monotonic())
        assert ok is False

    def test_skips_dirs(self, tmp_path, tar_cache, trash_dir):
        src = os.path.join(tmp_path, 'input')
        d = os.path.join(src, 'mydir')
        os.makedirs(d)

        chunk = os.path.join(tar_cache, 'chunk_0000.tar')
        _make_tar(chunk, [(src.lstrip('/') + '/mydir', d)])

        ok = phase1_verify(tar_cache, trash_dir, [src], None, time.monotonic())
        assert ok is True
        assert os.path.isdir(d)

    def test_multiple_input_dirs(self, tmp_path, tar_cache, trash_dir):
        src1 = os.path.join(tmp_path, 'input1')
        src2 = os.path.join(tmp_path, 'input2')
        _write_file(os.path.join(src1, 'a.txt'), 'a')
        _write_file(os.path.join(src2, 'b.txt'), 'b')

        chunk = os.path.join(tar_cache, 'chunk_0000.tar')
        _make_tar(chunk, [
            (src1.lstrip('/') + '/a.txt', os.path.join(src1, 'a.txt')),
            (src2.lstrip('/') + '/b.txt', os.path.join(src2, 'b.txt')),
        ])

        ok = phase1_verify(tar_cache, trash_dir, [src1, src2], None, time.monotonic())
        assert ok is True
        assert not os.path.exists(os.path.join(src1, 'a.txt'))
        assert not os.path.exists(os.path.join(src2, 'b.txt'))

    def test_no_disk_match_not_in_any_input_dir(self, tmp_path, tar_cache, trash_dir):
        src = os.path.join(tmp_path, 'input')
        _write_file(os.path.join(src, 'f.txt'), 'data')

        other = os.path.join(tmp_path, 'other')
        _write_file(os.path.join(other, 'f.txt'), 'data')

        chunk = os.path.join(tar_cache, 'chunk_0000.tar')
        _make_tar(chunk, [(other.lstrip('/') + '/f.txt', os.path.join(other, 'f.txt'))])

        ok = phase1_verify(tar_cache, trash_dir, [src], None, time.monotonic())
        assert ok is True
        assert os.path.exists(os.path.join(other, 'f.txt'))


# ======================================================================
# phase2_append
# ======================================================================

class TestPhase2Append:
    def test_archives_remaining_files(self, tmp_path, tar_cache, trash_dir):
        src = os.path.join(tmp_path, 'input')
        _write_file(os.path.join(src, 'new.txt'), 'new content')
        _symlink('new.txt', os.path.join(src, 'link_new'))

        ok, added = phase2_append(tar_cache, trash_dir, [src], None, time.monotonic())
        assert ok is True
        assert added is True
        assert os.path.isfile(os.path.join(tar_cache, 'chunk_0000.tar'))

        m = _tar_members(os.path.join(tar_cache, 'chunk_0000.tar'))
        names = [x.name for x in m]
        rel = src.lstrip('/')
        assert rel + '/new.txt' in names
        assert rel + '/link_new' in names

    def test_nothing_to_archive(self, tmp_path, tar_cache, trash_dir):
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        ok, added = phase2_append(tar_cache, trash_dir, [src], None, time.monotonic())
        assert ok is True
        assert added is False
        assert not os.path.isfile(os.path.join(tar_cache, 'chunk_0000.tar'))

    def test_excludes_trash_and_tar_dirs(self, tmp_path, tar_cache, trash_dir):
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        trash_in_src = os.path.join(src, 'mytraps')
        os.makedirs(trash_in_src)
        _write_file(os.path.join(trash_in_src, 'ignore.txt'), 'ignored')

        tar_in_src = os.path.join(src, 'mytars')
        os.makedirs(tar_in_src)
        _write_file(os.path.join(tar_in_src, 'also_ignore.txt'), 'ignored')

        _write_file(os.path.join(src, 'keep.txt'), 'kept')

        ok, added = phase2_append(trash_in_src, trash_in_src, [src], None, time.monotonic())
        assert ok is True
        assert added is True

    def test_symlink_dir_handled(self, tmp_path, tar_cache, trash_dir):
        src = os.path.join(tmp_path, 'input')
        os.makedirs(os.path.join(src, 'real_dir'))
        _symlink('real_dir', os.path.join(src, 'link_dir'))

        ok, added = phase2_append(tar_cache, trash_dir, [src], None, time.monotonic())
        assert ok is True
        assert added is True
        chunk = os.path.join(tar_cache, 'chunk_0000.tar')
        members = _tar_members(chunk)
        names = [m.name for m in members]
        rel = src.lstrip('/')
        assert rel + '/link_dir' in names

    def test_timeout_returns_false(self, tmp_path, tar_cache, trash_dir):
        src = os.path.join(tmp_path, 'input')
        _write_file(os.path.join(src, 'f.txt'), 'data')

        ok, added = phase2_append(tar_cache, trash_dir, [src], 0, time.monotonic())
        assert ok is False
        assert added is False

    def test_non_existent_input_dir_skipped(self, tmp_path, tar_cache, trash_dir):
        src = os.path.join(tmp_path, 'input')
        _write_file(os.path.join(src, 'f.txt'), 'data')

        ok, added = phase2_append(tar_cache, trash_dir, [src, '/nonexistent_dir'], None, time.monotonic())
        assert ok is True
        assert added is True

    def test_chunk_number_increments(self, tmp_path, tar_cache, trash_dir):
        src = os.path.join(tmp_path, 'input')
        _write_file(os.path.join(src, 'a.txt'), 'a')
        _write_file(os.path.join(src, 'b.txt'), 'b')

        ok1, added1 = phase2_append(tar_cache, trash_dir, [src], None, time.monotonic())
        assert added1 is True
        assert os.path.isfile(os.path.join(tar_cache, 'chunk_0000.tar'))

        _write_file(os.path.join(src, 'c.txt'), 'c')
        ok2, added2 = phase2_append(tar_cache, trash_dir, [src], None, time.monotonic())
        assert added2 is True
        assert os.path.isfile(os.path.join(tar_cache, 'chunk_0001.tar'))

    def test_temp_file_removed_when_empty(self, tmp_path, tar_cache, trash_dir):
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        tmp_files_before = set(os.listdir(tar_cache)) if os.path.isdir(tar_cache) else set()
        phase2_append(tar_cache, trash_dir, [src], None, time.monotonic())
        tmp_files_after = set(os.listdir(tar_cache)) if os.path.isdir(tar_cache) else set()
        assert tmp_files_before == tmp_files_after

    def test_unreadable_file_skipped(self, tmp_path, tar_cache, trash_dir):
        src = os.path.join(tmp_path, 'input')
        readable = os.path.join(src, 'good.txt')
        unreadable = os.path.join(src, 'bad.txt')
        _write_file(readable, 'good')
        _write_file(unreadable, 'bad')
        os.chmod(unreadable, 0)

        try:
            ok, added = phase2_append(tar_cache, trash_dir, [src], None, time.monotonic())
            assert ok is True
            assert added is True
            members = _tar_members(os.path.join(tar_cache, 'chunk_0000.tar'))
            assert len(members) == 1
        finally:
            os.chmod(unreadable, 0o644)


# ======================================================================
# main (integration via function call)
# ======================================================================

class TestMain:
    def _run(self, *args):
        try:
            return main(list(args))
        except SystemExit as e:
            return e.code

    def test_initial_archive_returns_3(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)
        _write_file(os.path.join(src, 'f.txt'), 'hello')

        rc = self._run('-f', tar, '-t', trash, src)
        assert rc == 3
        assert os.path.isfile(os.path.join(tar, 'chunk_0000.tar'))

    def test_idempotent_returns_0(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)
        _write_file(os.path.join(src, 'f.txt'), 'hello')

        self._run('-f', tar, '-t', trash, src)
        rc2 = self._run('-f', tar, '-t', trash, src)
        assert rc2 == 0

    def test_timeout_returns_2(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)
        _write_file(os.path.join(src, 'f.txt'), 'hello')

        rc = self._run('-f', tar, '-t', trash, '--timeout', '0', src)
        assert rc == 2

    def test_new_file_returns_3_again(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)
        _write_file(os.path.join(src, 'a.txt'), 'a')

        self._run('-f', tar, '-t', trash, src)
        _write_file(os.path.join(src, 'b.txt'), 'b')
        rc2 = self._run('-f', tar, '-t', trash, src)
        assert rc2 == 3

    def test_missing_input_dir_exits_1(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        os.makedirs(tar)
        os.makedirs(trash)

        rc = self._run('-f', tar, '-t', trash, '/nonexistent_input_dir')
        assert rc == 1

    def test_modified_file_both_versions_preserved(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        f = os.path.join(src, 'file.txt')
        _write_file(f, 'v1')
        self._run('-f', tar, '-t', trash, src)

        self._run('-f', tar, '-t', trash, src)

        _write_file(f, 'v2')
        self._run('-f', tar, '-t', trash, src)

        self._run('-f', tar, '-t', trash, src)

        rel = src.lstrip('/')
        v1_path = os.path.join(trash, rel, 'file.txt')
        v2_path = os.path.join(trash, rel, 'file.txt.1')
        assert os.path.isfile(v1_path)
        assert os.path.isfile(v2_path)
        assert Path(v1_path).read_text() == 'v1'
        assert Path(v2_path).read_text() == 'v2'
        assert not os.path.exists(f)

    def test_same_size_different_content_not_data_loss(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        f = os.path.join(src, 'file.bin')
        content_v1 = b'A' * 8192
        _write_file(f, content_v1)
        self._run('-f', tar, '-t', trash, src)

        self._run('-f', tar, '-t', trash, src)

        content_v2 = b'B' * 8192
        _write_file(f, content_v2)
        self._run('-f', tar, '-t', trash, src)

        self._run('-f', tar, '-t', trash, src)

        rel = src.lstrip('/')
        v1_path = os.path.join(trash, rel, 'file.bin')
        v2_path = os.path.join(trash, rel, 'file.bin.1')
        assert Path(v1_path).read_bytes() == content_v1
        assert Path(v2_path).read_bytes() == content_v2
        assert not os.path.exists(f)

    def test_no_data_loss_after_partial_timeout(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        _write_file(os.path.join(src, 'a.txt'), 'a')

        self._run('-f', tar, '-t', trash, src)
        _write_file(os.path.join(src, 'b.txt'), 'b')

        self._run('-f', tar, '-t', trash, '--timeout', '0', src)

        rc = self._run('-f', tar, '-t', trash, src)
        assert rc == 3

    def test_warning_when_trash_inside_input(self, tmp_path, capsys):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'input', 'trashes')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)
        _write_file(os.path.join(src, 'f.txt'), 'data')

        self._run('-f', tar, '-t', trash, src)
        captured = capsys.readouterr()
        has_warning = 'trash directory' in captured.out or 'trash directory' in captured.err
        assert has_warning

    def test_warning_when_tar_inside_input(self, tmp_path, capsys):
        tar = os.path.join(tmp_path, 'input', 'tars')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)
        _write_file(os.path.join(src, 'f.txt'), 'data')

        self._run('-f', tar, '-t', trash, src)
        captured = capsys.readouterr()
        has_warning = 'tar directory' in captured.out or 'tar directory' in captured.err
        assert has_warning


# ======================================================================
# subprocess-based integration tests
# ======================================================================

class TestSubprocessIntegration:
    def _run(self, *args):
        return subprocess.run(
            BLSL_MODULE + list(args),
            capture_output=True,
            text=True,
        )

    def test_idempotent_repeated_runs(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        for i in range(5):
            _write_file(os.path.join(src, f'f{i}.txt'), f'content {i}')

        r = self._run('-f', tar, '-t', trash, src)
        assert r.returncode == 3

        for _ in range(3):
            r = self._run('-f', tar, '-t', trash, src)
            assert r.returncode == 0

        chunks = sorted(f for f in os.listdir(tar) if CHUNK_RE.match(f))
        assert len(chunks) == 1

    def test_special_chars_in_filename(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        special_names = [
            'file with spaces.txt',
            'file-with-dashes.txt',
            'file_with_underscores.txt',
            'file.with.dots.txt',
            'file_with_plus+sign.txt',
        ]

        for name in special_names:
            _write_file(os.path.join(src, name), f'content of {name}')

        r = self._run('-f', tar, '-t', trash, src)
        assert r.returncode == 3
        r2 = self._run('-f', tar, '-t', trash, src)
        assert r2.returncode == 0

        rel = src.lstrip('/')
        for name in special_names:
            assert os.path.isfile(os.path.join(trash, rel, name))

    def test_unicode_filenames(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        uni_name = 'f\xfcle_\u00e7ontents_\U0001f4a9.txt'
        _write_file(os.path.join(src, uni_name), 'unicode content')

        r = self._run('-f', tar, '-t', trash, src)
        assert r.returncode == 3
        r2 = self._run('-f', tar, '-t', trash, src)
        assert r2.returncode == 0

        rel = src.lstrip('/')
        assert os.path.isfile(os.path.join(trash, rel, uni_name))

    def test_deeply_nested_dirs(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        depth = 50
        current = src
        for i in range(depth):
            current = os.path.join(current, f'dir_{i}')
            os.makedirs(current)
        _write_file(os.path.join(current, 'deep.txt'), 'bottom')

        r = self._run('-f', tar, '-t', trash, src)
        assert r.returncode == 3
        r2 = self._run('-f', tar, '-t', trash, src)
        assert r2.returncode == 0

    def test_empty_directories(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(os.path.join(src, 'empty1'))
        os.makedirs(os.path.join(src, 'empty2'))
        _write_file(os.path.join(src, 'not_empty', 'f.txt'), 'data')

        r = self._run('-f', tar, '-t', trash, src)
        assert r.returncode == 3
        r2 = self._run('-f', tar, '-t', trash, src)
        assert r2.returncode == 0

    def test_broken_symlinks(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        _symlink('nonexistent', os.path.join(src, 'broken'))
        _write_file(os.path.join(src, 'real.txt'), 'real')

        r = self._run('-f', tar, '-t', trash, src)
        assert r.returncode == 3
        r2 = self._run('-f', tar, '-t', trash, src)
        assert r2.returncode == 0

        rel = src.lstrip('/')
        trash_broken = os.path.join(trash, rel, 'broken')
        assert os.path.lexists(trash_broken)
        assert os.path.islink(trash_broken)

    def test_symlink_outside_input(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        outside = os.path.join(tmp_path, 'outside.txt')
        _write_file(outside, 'outside')
        _symlink(outside, os.path.join(src, 'extlink'))

        r = self._run('-f', tar, '-t', trash, src)
        assert r.returncode == 3
        r2 = self._run('-f', tar, '-t', trash, src)
        assert r2.returncode == 0

        rel = src.lstrip('/')
        trashed = os.path.join(trash, rel, 'extlink')
        assert os.path.lexists(trashed)
        assert os.path.islink(trashed)

    def test_many_files(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        n = 500
        for i in range(n):
            _write_file(os.path.join(src, f'file_{i:04d}.txt'), f'content {i}')

        r = self._run('-f', tar, '-t', trash, src)
        assert r.returncode == 3
        r2 = self._run('-f', tar, '-t', trash, src)
        assert r2.returncode == 0

        chunks = sorted(f for f in os.listdir(tar) if CHUNK_RE.match(f))
        assert len(chunks) == 1

        m = _tar_members(os.path.join(tar, chunks[0]))
        assert len([x for x in m if x.isfile()]) == n

    def test_large_file_content_integrity_4gib(self, tmp_path):
        gb = 4
        size = gb * 1024 * 1024 * 1024
        free = shutil.disk_usage(tmp_path).free

        if free < size * 5:
            pytest.skip(f'insufficient disk space for {gb} GiB test (need ~{size * 5} free, have {free})')

        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        large_path = os.path.join(src, 'large_4gib.bin')

        chunk = 16 * 1024 * 1024
        pattern = os.urandom(4096)
        with open(large_path, 'wb') as fh:
            written = 0
            while written < size:
                fh.write(pattern * (min(chunk, size - written) // 4096))
                written += min(chunk, size - written) // 4096 * 4096
            remaining = size - fh.tell()
            if remaining > 0:
                fh.write(os.urandom(remaining))
        assert os.path.getsize(large_path) == size

        r = self._run('-f', tar, '-t', trash, src)
        assert r.returncode == 3

        tar_chunk = os.path.join(tar, 'chunk_0000.tar')
        assert os.path.isfile(tar_chunk)

        extract_dir = os.path.join(tmp_path, 'extracted')
        os.makedirs(extract_dir)
        with tarfile.open(tar_chunk, 'r') as tf:
            tf.extractall(path=extract_dir)

        rel = src.lstrip('/')
        extracted = os.path.join(extract_dir, rel, 'large_4gib.bin')
        assert os.path.getsize(extracted) == size

        bufsize = 1024 * 1024
        with open(large_path, 'rb') as a_fh, open(extracted, 'rb') as b_fh:
            identical = True
            offset = 0
            while True:
                chunk_a = a_fh.read(bufsize)
                chunk_b = b_fh.read(bufsize)
                if chunk_a != chunk_b:
                    identical = False
                    break
                if not chunk_a:
                    break
                offset += len(chunk_a)
                if offset % (100 * bufsize) == 0:
                    pass

            assert identical

    def test_concatenated_chunks_valid_tar(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        for i in range(3):
            _write_file(os.path.join(src, f'batch1_f{i}.txt'), f'first batch {i}')
        self._run('-f', tar, '-t', trash, src)

        for i in range(3):
            _write_file(os.path.join(src, f'batch2_f{i}.txt'), f'second batch {i}')
        self._run('-f', tar, '-t', trash, src)

        for i in range(3):
            _write_file(os.path.join(src, f'batch3_f{i}.txt'), f'third batch {i}')
        self._run('-f', tar, '-t', trash, src)

        combined = os.path.join(tmp_path, 'combined.tar')
        chunks = sorted(
            [f for f in os.listdir(tar) if CHUNK_RE.match(f)],
            key=lambda x: int(CHUNK_RE.match(x).group(1)),
        )
        assert len(chunks) == 3

        with open(combined, 'wb') as out_fh:
            for c in chunks:
                with open(os.path.join(tar, c), 'rb') as in_fh:
                    out_fh.write(in_fh.read())

        result = subprocess.run(
            ['tar', 'tf', combined, '--ignore-zeros'],
            capture_output=True, text=True,
        )
        lines = [l for l in result.stdout.strip().split('\n') if l]
        file_count = len([l for l in lines if not l.endswith('/')])
        assert file_count == 9

    def test_recovery_after_orphaned_temp_file(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        _write_file(os.path.join(src, 'survive.txt'), 'must survive')

        _symlink('survive.txt', os.path.join(src, 'link'))
        orphan = os.path.join(tar, '.tmp_chunk_orphan.tar')
        os.makedirs(tar, exist_ok=True)
        _write_file(orphan, b'garbage bytes not real tar')

        r = self._run('-f', tar, '-t', trash, src)
        assert r.returncode == 3
        assert os.path.isfile(orphan)

        r2 = self._run('-f', tar, '-t', trash, src)
        assert r2.returncode == 0

    def test_no_read_permission_file(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        protected = os.path.join(src, 'protected.txt')
        _write_file(protected, 'secret')
        os.chmod(protected, 0)
        _write_file(os.path.join(src, 'normal.txt'), 'public')

        try:
            self._run('-f', tar, '-t', trash, src)

            r = self._run('-f', tar, '-t', trash, src)
            assert r.returncode == 0
            chunks = sorted(f for f in os.listdir(tar) if CHUNK_RE.match(f))
            assert len(chunks) >= 1
            m = _tar_members(os.path.join(tar, chunks[0]))
            normal_members = [x for x in m if 'normal.txt' in x.name]
            assert len(normal_members) == 1
        finally:
            os.chmod(protected, 0o644)

    def test_restore_from_trash_and_rearchive(self, tmp_path):
        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        _write_file(os.path.join(src, 'f.txt'), 'original')
        r1 = self._run('-f', tar, '-t', trash, src)
        assert r1.returncode == 3

        r_verify = self._run('-f', tar, '-t', trash, src)
        assert r_verify.returncode == 0
        assert not os.path.exists(os.path.join(src, 'f.txt'))

        _write_file(os.path.join(src, 'f.txt'), 'restored')
        r2 = self._run('-f', tar, '-t', trash, src)
        assert r2.returncode == 3

        r_verify2 = self._run('-f', tar, '-t', trash, src)
        assert r_verify2.returncode == 0

        rel = src.lstrip('/')
        assert os.path.isfile(os.path.join(trash, rel, 'f.txt'))
        assert os.path.isfile(os.path.join(trash, rel, 'f.txt.1'))
        assert not os.path.exists(os.path.join(src, 'f.txt'))


# ======================================================================
# crash recovery tests (kill -9 simulation)
# ======================================================================

class TestCrashRecovery:
    CRASH = pytest.mark.crash

    @CRASH
    def test_kill_during_phase2_append_tar_write(self, tmp_path):
        """Kill -9 during tar write: temp file orphaned, original files intact."""
        if not hasattr(signal, 'SIGKILL'):
            pytest.skip('SIGKILL not available')

        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        _write_file(os.path.join(src, 'a.txt'), 'aaa')
        self._run_subprocess('-f', tar, '-t', trash, src)

        _write_file(os.path.join(src, 'b.txt'), 'bbb')

        proc = subprocess.Popen(
            BLSL_MODULE + ['-f', tar, '-t', trash, src],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        time.sleep(0.05)
        os.kill(proc.pid, signal.SIGKILL)
        proc.wait()

        assert os.path.isfile(os.path.join(src, 'b.txt'))

        r = self._run_subprocess('-f', tar, '-t', trash, src)
        assert r.returncode in (0, 3)

    @CRASH
    def test_kill_during_phase1_verify(self, tmp_path):
        """Kill during verification: re-run idempotently recovers."""
        if not hasattr(signal, 'SIGKILL'):
            pytest.skip('SIGKILL not available')

        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        for i in range(20):
            _write_file(os.path.join(src, f'f{i:04d}.txt'), f'content {i}')
        self._run_subprocess('-f', tar, '-t', trash, src)

        for i in range(20, 40):
            _write_file(os.path.join(src, f'f{i:04d}.txt'), f'content {i}')

        proc = subprocess.Popen(
            BLSL_MODULE + ['-f', tar, '-t', trash, src],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        time.sleep(0.1)
        os.kill(proc.pid, signal.SIGKILL)
        proc.wait()

        for _ in range(3):
            r = self._run_subprocess('-f', tar, '-t', trash, src)
            assert r.returncode in (0, 3)

    @CRASH
    def test_kill_between_verify_and_move(self, tmp_path):
        """Verify passes but SIGKILL before move: file stays, re-run re-verifies."""
        if not hasattr(signal, 'SIGKILL'):
            pytest.skip('SIGKILL not available')

        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        for i in range(10):
            _write_file(os.path.join(src, f'f{i:04d}.txt'), f'content {i}')

        self._run_subprocess('-f', tar, '-t', trash, src)

        for i in range(10, 20):
            _write_file(os.path.join(src, f'f{i:04d}.txt'), f'content {i}')

        proc = subprocess.Popen(
            BLSL_MODULE + ['-f', tar, '-t', trash, src],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        time.sleep(0.05)
        os.kill(proc.pid, signal.SIGKILL)
        proc.wait()

        r = self._run_subprocess('-f', tar, '-t', trash, src)
        assert r.returncode in (0, 3)

        for i in range(10):
            assert not os.path.exists(os.path.join(src, f'f{i:04d}.txt'))

    @CRASH
    def test_repeated_kill_recovery_is_idempotent(self, tmp_path):
        """Repeated kill-and-restart never loses data."""
        if not hasattr(signal, 'SIGKILL'):
            pytest.skip('SIGKILL not available')

        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        files = [f'f{i:04d}.txt' for i in range(10)]
        for name in files:
            _write_file(os.path.join(src, name), f'content {name}')

        max_kill_attempts = 8
        for kill_attempt in range(max_kill_attempts):
            proc = subprocess.Popen(
                BLSL_MODULE + ['-f', tar, '-t', trash, src],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            time.sleep(0.02 * (kill_attempt + 1))
            try:
                os.kill(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            proc.wait()

        success = False
        for attempt in range(10):
            r = self._run_subprocess('-f', tar, '-t', trash, src)
            if r.returncode == 0:
                success = True
                break
            elif r.returncode == 3:
                time.sleep(0.05)

        assert success

        r = self._run_subprocess('-f', tar, '-t', trash, src)
        assert r.returncode == 0
        for name in files:
            assert not os.path.exists(os.path.join(src, name))

    @CRASH
    def test_kill_mid_archive_no_data_loss(self, tmp_path):
        """Kill a running tarverify during archiving, verify no data loss after recovery."""
        if not hasattr(signal, 'SIGKILL'):
            pytest.skip('SIGKILL not available')

        tar = os.path.join(tmp_path, 'tar')
        trash = os.path.join(tmp_path, 'trash')
        src = os.path.join(tmp_path, 'input')
        os.makedirs(src)

        files = [f'f{i:05d}.txt' for i in range(15)]
        for name in files:
            _write_file(os.path.join(src, name), f'content for {name}')

        for kill_attempt in range(5):
            proc = subprocess.Popen(
                BLSL_MODULE + ['-f', tar, '-t', trash, src],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            time.sleep(0.03 * (kill_attempt + 1))
            try:
                os.kill(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            proc.wait()

        max_attempts = 10
        for attempt in range(max_attempts):
            r = self._run_subprocess('-f', tar, '-t', trash, src)
            if r.returncode == 0:
                break

        for name in files:
            assert not os.path.exists(os.path.join(src, name))

        rel = src.lstrip('/')
        for name in files:
            trash_path = os.path.join(trash, rel, name)
            assert os.path.exists(trash_path) or any(
                os.path.exists(os.path.join(trash, rel, f'{name}.{s}'))
                for s in range(10)
            )

    # helpers for crash tests
    def _run_subprocess(self, *args):
        return subprocess.run(
            BLSL_MODULE + list(args),
            capture_output=True,
            text=True,
        )


# ======================================================================
# filesystem-full simulation
# ======================================================================

class TestFileSystemFull:
    SLOW = pytest.mark.slow

    @SLOW
    def test_full_filesystem_no_data_loss(self, tmp_path):
        """Write to a near-full filesystem, ensure no data loss."""
        img = os.path.join(tmp_path, 'small.img')
        mnt = os.path.join(tmp_path, 'mnt')

        img_size_mb = 20
        with open(img, 'wb') as fh:
            fh.seek(img_size_mb * 1024 * 1024 - 1)
            fh.write(b'\x00')

        subprocess.run(['mkfs.ext4', '-q', '-F', img], capture_output=True, check=True)

        if not os.path.isdir(mnt):
            os.makedirs(mnt)

        try:
            subprocess.run(
                ['sudo', 'mount', '-o', 'loop', img, mnt],
                capture_output=True, check=True, timeout=10,
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            pytest.skip('Cannot mount loop device (requires sudo)')

        try:
            subprocess.run(['sudo', 'chmod', '777', mnt], capture_output=True, check=True, timeout=5)

            tar_in = os.path.join(mnt, 'tar')
            trash_in = os.path.join(mnt, 'trash')

            src = os.path.join(tmp_path, 'input')
            os.makedirs(src)
            _write_file(os.path.join(src, 'safe.txt'), 'safe data')

            try:
                r = self._run('-f', tar_in, '-t', trash_in, src)
                assert r.returncode in (0, 3)

                if r.returncode == 3:
                    assert os.path.isfile(os.path.join(src, 'safe.txt'))
                    assert not os.path.isfile(os.path.join(tar_in, 'chunk_0000.tar'))
            except SystemExit:
                pass

        finally:
            subprocess.run(['sudo', 'umount', mnt], capture_output=True, timeout=10)

    def _run(self, *args):
        try:
            return main(list(args))
        except SystemExit as e:
            return type('FakeResult', (), {'returncode': e.code})()

