#!/usr/bin/env bash
set -euo pipefail

TARVERIFY="python3 -m blsl tarverify"
TDIR=$(mktemp -d)
trap 'rm -rf "$TDIR"' EXIT

# Short aliases for test paths
IDIR="$TDIR/input"
TRASH="$TDIR/trash"
TAR="$TDIR/tar"

echo "=== test_tarverify.sh ==="

# --- setup ---
mkdir -p "$IDIR/a" "$IDIR/b" "$TAR" "$TRASH"

echo "file_a" > "$IDIR/a/f1.txt"
echo "file_b" > "$IDIR/b/f2.txt"
ln -s f1.txt "$IDIR/a/link1.txt"

# ===== test 1: first run — archive everything =====
echo "--- Test 1: initial archive ---"
$TARVERIFY -f "$TAR" -t "$TRASH" "$IDIR/a" "$IDIR/b"

# Absolute-encoded paths under trash
[ -f "$TRASH$IDIR/a/f1.txt" ] || { echo "FAIL: f1 not trashed at $TRASH$IDIR/a/f1.txt"; exit 1; }
[ -f "$TRASH$IDIR/b/f2.txt" ] || { echo "FAIL: f2 not trashed at $TRASH$IDIR/b/f2.txt"; exit 1; }
[ -L "$TRASH$IDIR/a/link1.txt" ] || { echo "FAIL: link1 not trashed at $TRASH$IDIR/a/link1.txt"; exit 1; }

CHUNKS=$(ls "$TAR"/chunk_*.tar 2>/dev/null | wc -l)
[ "$CHUNKS" -eq 1 ] || { echo "FAIL: expected 1 chunk, got $CHUNKS"; exit 1; }

echo "PASS"

# ===== test 2: re-run idempotent =====
echo "--- Test 2: re-run idempotent ---"
$TARVERIFY -f "$TAR" -t "$TRASH" "$IDIR/a" "$IDIR/b"

CHUNKS2=$(ls "$TAR"/chunk_*.tar 2>/dev/null | wc -l)
[ "$CHUNKS2" -eq 1 ] || { echo "FAIL: expected 1 chunk, got $CHUNKS2"; exit 1; }

echo "PASS"

# ===== test 3: new file after archive =====
echo "--- Test 3: new file after archive ---"
echo "newfile" > "$IDIR/a/new.txt"
$TARVERIFY -f "$TAR" -t "$TRASH" "$IDIR/a" "$IDIR/b"

[ -f "$TRASH$IDIR/a/new.txt" ] || { echo "FAIL: new.txt not trashed"; exit 1; }

CHUNKS3=$(ls "$TAR"/chunk_*.tar | wc -l)
[ "$CHUNKS3" -eq 2 ] || { echo "FAIL: expected 2 chunks, got $CHUNKS3"; exit 1; }

echo "PASS"

# ===== test 4: changed content — re-archived =====
echo "--- Test 4: changed content ---"
mkdir -p "$IDIR/c"
echo "v1" > "$IDIR/c/change.txt"
$TARVERIFY -f "$TAR" -t "$TRASH" "$IDIR/c"
[ -f "$TRASH$IDIR/c/change.txt" ] || { echo "FAIL: v1 not trashed"; exit 1; }

echo "v2" > "$IDIR/c/change.txt"
$TARVERIFY -f "$TAR" -t "$TRASH" "$IDIR/c"
[ -f "$TRASH$IDIR/c/change.txt.1" ] || { echo "FAIL: v2 not trashed as .1"; exit 1; }

echo "PASS"

# ===== test 5: timeout exit code 2 =====
echo "--- Test 5: timeout ---"
echo "slowfile" > "$IDIR/a/slow1.txt"
set +e
$TARVERIFY -f "$TAR" -t "$TRASH" --timeout 0 "$IDIR/a"
ret=$?
set -e
[ "$ret" -eq 2 ] || { echo "FAIL: expected exit 2, got $ret"; exit 1; }
echo "PASS"

# ===== test 6: concatenated chunks =====
echo "--- Test 6: cat chunks ---"
cat "$TAR"/chunk_*.tar > "$TDIR/combined.tar"
tar tf "$TDIR/combined.tar" > /dev/null && echo "combined tar OK"
echo "PASS"

echo "=== ALL TESTS PASSED ==="
