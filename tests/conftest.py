import os
import shutil
import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def tmp_dir():
    path = tempfile.mkdtemp(prefix='blsl_test_')
    yield path
    shutil.rmtree(path, ignore_errors=True)


@pytest.fixture
def input_tree(tmp_dir):
    src = os.path.join(tmp_dir, 'input')
    os.makedirs(src)
    os.makedirs(os.path.join(src, 'subdir', 'deep'))

    (Path(src) / 'file_a.txt').write_text('hello world\n')
    (Path(src) / 'empty_file.dat').write_text('')
    (Path(src) / 'binary.bin').write_bytes(b'\x00\x01\x02\x03\xff\xfe\xfd')

    big = b'A' * 100000 + b'B' * 100000
    (Path(src) / 'large_file.bin').write_bytes(big)

    deep_dir = os.path.join(src, 'subdir', 'deep')
    (Path(deep_dir) / 'deep_file.txt').write_text('deep content')

    os.makedirs(os.path.join(src, 'empty_dir'))
    os.symlink('file_a.txt', os.path.join(src, 'link_to_file'))
    os.symlink('subdir', os.path.join(src, 'link_to_dir'))
    os.symlink('nonexistent', os.path.join(src, 'broken_link'))

    yield src
