#!/usr/bin/env bash
# Golden tests: diff mytools against real bedtools.
# Usage: ./tests/run_golden.sh
set -uo pipefail

MYTOOLS=${MYTOOLS:-./bin/mytools}   # override to test a different build
DATA=$(dirname "$0")/../data
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
pass=0; fail=0

# report <name> <got_rc> <want_rc>
#   compares exit codes and $tmp/got against $tmp/want. Shared by check and
#   check_stdin so the two cannot drift apart.
report() {
  local name=$1 got_rc=$2 want_rc=$3

  if [[ $got_rc -ne $want_rc ]]; then
    echo "FAIL $name (exit $got_rc, bedtools gave $want_rc)"
    sed 's/^/      /' "$tmp/got.err" | head -3
    (( fail++ )); return
  fi
  if diff -q "$tmp/want" "$tmp/got" >/dev/null; then
    echo "ok   $name"; (( pass++ ))
  else
    echo "FAIL $name"
    diff -u "$tmp/want" "$tmp/got" | sed 's/^/      /' | head -20
    (( fail++ ))
  fi
}

# check <name> -- <args...>
#   runs "$MYTOOLS <args>" and "bedtools <args>", diffs them
check() {
  local name=$1; shift; shift        # drop the literal --
  "$MYTOOLS" "$@" > "$tmp/got"  2>"$tmp/got.err"
  local got_rc=$?
  bedtools   "$@" > "$tmp/want" 2>/dev/null
  local want_rc=$?
  report "$name" $got_rc $want_rc
}

# check_stdin <name> <file> -- <args...>
#   same, but <file> is piped in on stdin -- for the cases where an input is `-`.
check_stdin() {
  local name=$1 input=$2; shift 2; shift
  "$MYTOOLS" "$@" < "$input" > "$tmp/got"  2>"$tmp/got.err"
  local got_rc=$?
  bedtools   "$@" < "$input" > "$tmp/want" 2>/dev/null
  local want_rc=$?
  report "$name" $got_rc $want_rc
}

# check_rc <name> <expected_rc> -- <args...>
#   asserts an exit code without consulting bedtools, and that stdout stayed empty.
#   For usage errors only: SPEC.md section 7 has us exit 2 where bedtools exits 1,
#   so the oracle does not adjudicate these and there is nothing to diff.
check_rc() {
  local name=$1 want_rc=$2; shift 2; shift
  "$MYTOOLS" "$@" > "$tmp/got" 2>"$tmp/got.err"
  local got_rc=$?

  if [[ $got_rc -ne $want_rc ]]; then
    echo "FAIL $name (exit $got_rc, expected $want_rc)"
    sed 's/^/      /' "$tmp/got.err" | head -3
    (( fail++ )); return
  fi
  if [[ -s "$tmp/got" ]]; then
    echo "FAIL $name (exit $got_rc as expected, but wrote to stdout)"
    sed 's/^/      /' "$tmp/got" | head -3
    (( fail++ )); return
  fi
  echo "ok   $name"; (( pass++ ))
}

# --- sort (#5) -- SPEC.md section 8 case 1 -------------------------------------
# a.bed is deliberately unsorted, and the answer includes the zero-length features
# a07 (chr1 500 500) and a12 (chr2 0 0). Whatever bedtools prints is correct.
check "sort a.bed" -- sort -i "$DATA/a.bed"

# Same case through stdin. `-` must behave exactly like the file argument.
check_stdin "sort a.bed (stdin)" "$DATA/a.bed" -- sort -i -

# Exit code on bad data, diffed against the oracle: both must print nothing to
# stdout and exit 1 (SPEC.md section 7). The message text is our own; stderr is
# not part of the contract.
printf 'chr1\tfoo\t200\n' > "$tmp/malformed.bed"
check "sort malformed input" -- sort -i "$tmp/malformed.bed"

# --- intersect (#7) -- SPEC.md section 8 cases 5-10 -----------------------------
# No sorting needed: intersect indexes -b and streams -a.
check "intersect"     -- intersect     -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -u"  -- intersect -u  -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -v"  -- intersect -v  -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -wa" -- intersect -wa -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -wb" -- intersect -wb -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -s"  -- intersect -s  -a "$DATA/a.bed" -b "$DATA/b.bed"

# The same answer as case 5, reached with -a on stdin.
check_stdin "intersect -a - (stdin)" "$DATA/a.bed" -- intersect -a - -b "$DATA/b.bed"

# -u, -v, -wa and -wb are mutually exclusive (SPEC.md section 5). bedtools accepts
# some of these combinations, so these are ours to assert, not the oracle's.
check_rc "intersect -u -v rejected"   2 -- intersect -u -v   -a "$DATA/a.bed" -b "$DATA/b.bed"
check_rc "intersect -u -wa rejected"  2 -- intersect -u -wa  -a "$DATA/a.bed" -b "$DATA/b.bed"
check_rc "intersect -wa -wb rejected" 2 -- intersect -wa -wb -a "$DATA/a.bed" -b "$DATA/b.bed"

# Add the rest here. merge and closest need sorted input -- sort into $tmp first.
# The 17 cases are enumerated in SPEC.md section 8, and each subcommand issue
# (#5-#9) names the ones it owns.

echo "---"
echo "$pass passed, $fail failed"
[[ $fail -eq 0 ]]
