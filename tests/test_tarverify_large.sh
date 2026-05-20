#!/usr/bin/env bash
# Large-file test for tarverify.py.
# Set TARVERIFY_LARGE_GB=10 for the full 10 GiB test (default: 10).
set -euo pipefail

TARVERIFY="python3 -m blsl tarverify"

SIZE_GB="${TARVERIFY_LARGE_GB:-10}"
SIZE_BYTES=$((SIZE_GB * 1024 * 1024 * 1024))
SIZE_MB=$((SIZE_GB * 1024))

mkdir -p tmp
TDIR=$(mktemp -d -p $PWD/tmp)
trap 'rm -rf "$TDIR"' EXIT

IDIR="$TDIR/input"
TRASH="$TDIR/trash"
TAR="$TDIR/tar"
mkdir -p "$IDIR" "$TRASH" "$TAR"

# ---- space check ----
FREE_KB=$(df --output=avail "$TDIR" 2>/dev/null | tail -1)
NEED_KB=$((SIZE_GB * 1024 * 1024 * 3))
if [ "$FREE_KB" -lt "$NEED_KB" ]; then
    echo "SKIP: need ~${NEED_KB} KB free, have ${FREE_KB} KB"
    exit 0
fi

echo "=== test_tarverify_large.sh (${SIZE_GB} GiB files) ==="

# ---- helper: create a large file with known content ----
create_large() {
    local path="$1" pattern="$2"
    # Write a header we can easily verify, then fill the rest with the pattern
    printf "MAGIC_%s\n" "$pattern" > "$path"
    (yes "$pattern" || true) | dd bs=1M count=$SIZE_MB iflag=fullblock >> "$path" 2>/dev/null
}

echo "--- Creating ${SIZE_GB} GiB file (this may take a minute) ---"
create_large "$IDIR/large_a.dat" "AAAA"
echo "--- Creating second ${SIZE_GB} GiB file ---"
create_large "$IDIR/large_b.dat" "BBBB"

echo "--- Test 1: initial archive of large files ---"
$TARVERIFY -f "$TAR" -t "$TRASH" "$IDIR"

# Both should be trashed
[ -f "$TRASH$IDIR/large_a.dat" ] || { echo "FAIL: large_a.dat not trashed"; exit 1; }
[ -f "$TRASH$IDIR/large_b.dat" ] || { echo "FAIL: large_b.dat not trashed"; exit 1; }
echo "PASS: both files archived and trashed"

echo "--- Test 2: re-run idempotent with large files ---"
$TARVERIFY -f "$TAR" -t "$TRASH" "$IDIR"
CHUNKS=$(ls "$TAR"/chunk_*.tar | wc -l)
[ "$CHUNKS" -eq 1 ] || { echo "FAIL: expected 1 chunk after idempotent re-run, got $CHUNKS"; exit 1; }
echo "PASS: no new chunk created"

echo "--- Test 3: changed large file triggers re-archive ---"
# New version on disk while old version remains in trash → .1 suffix
create_large "$IDIR/large_a.dat" "ZZZZ"

$TARVERIFY -f "$TAR" -t "$TRASH" "$IDIR"

# Should be re-archived in a new chunk with .1 suffix (trash clash)
[ -f "$TRASH$IDIR/large_a.dat.1" ] || { echo "FAIL: changed large_a.dat not re-archived as .1"; exit 1; }
CHUNKS2=$(ls "$TAR"/chunk_*.tar | wc -l)
[ "$CHUNKS2" -eq 2 ] || { echo "FAIL: expected 2 chunks, got $CHUNKS2"; exit 1; }
echo "PASS: changed large file re-archived"

echo "--- Test 4: verify concatenated chunks ---"
cat "$TAR"/chunk_*.tar > "$TDIR/combined.tar"
exp_size=0; for f in "$TAR"/chunk_*.tar; do exp_size=$((exp_size + $(stat -c%s "$f"))); done
act_size=$(stat -c%s "$TDIR/combined.tar")
[ "$exp_size" -eq "$act_size" ] || { echo "FAIL: combined size mismatch"; exit 1; }
member_count=$(tar tf "$TDIR/combined.tar" --ignore-zeros 2>/dev/null | wc -l)
[ "$member_count" -ge 3 ] || { echo "FAIL: expected 3+ members in combined tar, got $member_count"; exit 1; }
echo "PASS: combined tar ($act_size bytes, $member_count members)"

echo "=== ALL LARGE FILE TESTS PASSED ==="
