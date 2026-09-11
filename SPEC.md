# SPEC.md — mytools

A small reimplementation of a subset of bedtools. Real `bedtools` is the oracle: if our
output differs from it on the same input, we are wrong.

Decisions here were made in the 0:35 brainstorm. Where the answer was "do it like
bedtools", the spec records what bedtools actually does — measured on this VM
(v2.31.1), not assumed.

---

## 1. Scope

v1 ships five subcommands: `sort`, `merge`, `intersect`, `subtract`, `closest`.

Explicitly **not** in v1:

- `-f` / `-F` fractional overlap, `-r`, `-e`
- `-c` counting, `-wo` / `-loj` output modes
- `merge -c` / `-o` column operations
- `sort -g` / `-faidx` genome-ordered sorting
- `closest -t` tie-breaking (default "all ties" behaviour only), `-k`, `-D`
- GFF3 and VCF input
- `-header`
- Multiple `-b` files

Scope is a decision. Five subcommands finished beats eight half-built ones.

## 2. Invocation

```
mytools sort      -i <file|->
mytools merge     -i <file|-> [-d N] [-s]
mytools intersect -a <file|-> -b <file> [-u|-v|-wa|-wb] [-s]
mytools subtract  -a <file|-> -b <file> [-A] [-s]
mytools closest   -a <file|-> -b <file> [-d] [-io] [-s]
mytools --version
```

Flag names and meanings match bedtools exactly. `-` means stdin; at most one input per
invocation may be `-`. `-b` must be a seekable file (we read it fully; see §6).

**Entry point.** `./mytools` at the repo root is an executable wrapper that imports the
package. Golden tests invoke `"${MYTOOLS:-./mytools}"` so the binary can be relocated
without editing every test.

## 3. Input formats

BED3 through BED12, tab-separated:

```
chrom  start  end  [name  score  strand  thickStart  thickEnd  itemRgb  blockCount  blockSizes  blockStarts]
```

- Column count may vary between lines. Trailing columns are preserved, not parsed —
  only `chrom`, `start`, `end` and `strand` (column 6) carry meaning in v1.
- Input may come from a file argument or stdin.
- Lines beginning with `#`, `track` or `browser`, and blank lines, are **skipped
  silently** — matching bedtools. They do not appear in output.
- Gzip input is supported and detected by the magic bytes `1f 8b` at the start of the
  stream, so it works for both files and stdin regardless of filename.
- `start` and `end` are non-negative integers. `start > end` is an error (§7).
- `start == end` (zero-length) is **legal**. See §4.

## 4. Interval semantics

BED is **0-based, half-open**. `chr1 100 200` covers bases 100..199.

- **Overlap predicate:** `a.start < b.end AND b.start < a.end`. Strict `<` on both
  sides. Every off-by-one bug in this project lives in that one line.
- **Minimum overlap:** 1bp. Any shared base counts; there is no fractional threshold in
  v1.
- **Bookended** intervals (`a.end == b.start`) do **not** overlap — the predicate above
  is false for them.
- Bookended intervals **do** merge at `-d 0` (the default), because `merge` joins
  features whose gap is `<= d`, and the gap between bookended features is 0.
- **Zero-length** intervals (`start == end`) are legal and appear in the fixtures
  (`data/b.bed` line 3 is `chr1 100 100`). Real bedtools handles them in ways nobody
  predicts from first principles. **Do not reason about them — run bedtools and match
  whatever it prints.** Encode the result with a comment saying it is oracle behaviour.

### Strand

`-s` means "same strand only" and applies to `merge`, `intersect`, `subtract` and
`closest`. Features with `.` or a missing strand column are treated as bedtools treats
them — run it and match. `data/a.bed` contains `a03`/`a04`, an identical span on
opposite strands, specifically to exercise this.

## 5. Output

- Tab-separated, LF line endings, trailing newline on the final line.
- Empty result: print nothing, exit 0.
- **`sort`** — all input columns preserved. Order is chrom (lexicographic, so `chr17`
  sorts before `chr7`), then `start`, then `end`.
- **`merge`** — BED3 only (`chrom start end`). Input columns are dropped. With `-s`,
  bedtools appends the strand column; match it.
- **`intersect`** — default prints the intersected region carrying `-a`'s trailing
  columns. `-wa` prints `-a`'s original interval, once per overlapping `-b` feature.
  `-wb` appends the overlapping `-b` feature. `-u` prints each `-a` feature at most
  once. `-v` prints `-a` features with no overlap. `-u`, `-v`, `-wa` and `-wb` are
  mutually exclusive.
- **`subtract`** — `-a` features with `-b` regions removed. One feature may become two,
  or vanish. `-A` removes the entire `-a` feature if it overlaps `-b` at all.
- **`closest`** — for each `-a` feature, the nearest `-b` feature(s), ties all reported.
  `-d` appends the distance as a final column (0 if they overlap). `-io` ignores
  overlapping features and reports the nearest non-overlapping one.
- Input order is preserved for `intersect`, `subtract` and `closest`. bedtools does not
  reorder `-a`, and neither do we.

## 6. Memory model

- **`sort`** holds the whole input in memory. It fundamentally cannot stream.
- **`merge`** streams and requires pre-sorted input. It does **not** sort for you;
  unsorted input is an error (§7), same as bedtools.
- **`closest`** requires both `-a` and `-b` pre-sorted, and errors otherwise.
- **`intersect`** and **`subtract`** load `-b` into memory (grouped by chrom, sorted by
  start, binary-searched) and stream `-a`.
- Target: inputs up to ~10⁶ intervals. That covers the ~500,000 intervals produced by
  `bedtools bamtobed` on the stretch-goal BAM.
- No mmap, no index files, no threads.

## 7. Errors and exit codes

Errors go to **stderr**. stdout carries data only, because it gets piped.

| Situation | exit | matches bedtools? |
|---|---|---|
| Success, including empty output | 0 | yes |
| Malformed line / non-integer coordinates | 1 | yes |
| `start > end` | 1 | yes |
| Unsorted input to `merge` or `closest` | 1 | yes |
| Unknown flag / missing required argument | 2 | **no** — bedtools exits 1 |
| Input file does not exist | 2 | **no** — bedtools exits 1 |
| `--version`, `--help` | 0 | yes |

**Why the split.** Measured on this VM, bedtools exits `1` for every error without
exception. We match it on everything a golden test compares — the data errors — so
those tests diff exit codes against the oracle for free. Usage errors get `2` because
no meaningful golden comparison exists for them (the message text differs anyway) and
the distinction is worth having. This is the only deliberate behavioural deviation in
v1, and it is deliberately confined to cases the oracle never adjudicates.

Error messages name the file and line where known:
`a.bed:14: start > end (500 > 400)`. Message text is not compared against bedtools.

## 8. Correctness

Real `bedtools` is the oracle. Non-negotiable.

**17 golden cases** — one per subcommand default, one per flag:

| # | Case |
|---|---|
| 1 | `sort -i data/a.bed` |
| 2 | `merge -i <sorted a.bed>` |
| 3 | `merge -d 10 -i <sorted a.bed>` |
| 4 | `merge -s -i <sorted a.bed>` |
| 5 | `intersect -a data/a.bed -b data/b.bed` |
| 6 | `intersect -u -a … -b …` |
| 7 | `intersect -v -a … -b …` |
| 8 | `intersect -wa -a … -b …` |
| 9 | `intersect -wb -a … -b …` |
| 10 | `intersect -s -a … -b …` |
| 11 | `subtract -a data/a.bed -b data/b.bed` |
| 12 | `subtract -A -a … -b …` |
| 13 | `subtract -s -a … -b …` |
| 14 | `closest -a <sorted a.bed> -b <sorted b.bed>` |
| 15 | `closest -d -a … -b …` |
| 16 | `closest -io -a … -b …` |
| 17 | `closest -s -a … -b …` |

Each case compares **stdout and exit code**. At least one case per subcommand must also
be exercised reading from **stdin** via `-`.

`merge` and `closest` cases sort into a temp file first — `data/a.bed` and `data/b.bed`
are deliberately unsorted.

Also required: `mytools --version` prints a version and exits 0; `mytools` with no
arguments prints usage to stderr and exits 2.

**Accepted deviations from bedtools: none**, other than the usage-error exit codes in
§7. If you find one you cannot fix, write it down here with the reason.

## 9. Language and layout

Python — every subcommand and every test. No second language without asking.

```
mytools                  executable wrapper, imports the package
mytools/
  __init__.py
  cli.py                 argument parsing, dispatch, exit codes
  bed.py                 parsing, the overlap predicate, shared interval logic
  sort.py
  merge.py
  intersect.py
  subtract.py
  closest.py
tests/
  run_golden.sh          the 17 cases above, diffed against bedtools
  test_*.py              pytest unit tests, one per edge case, no bedtools needed
```

- A module per subcommand, so five parallel agents don't collide in one file. Shared
  interval logic lives in `bed.py` and is the one file that needs coordination.
- Unit tests run under **pytest**, without bedtools installed, in seconds.
- No third-party runtime dependencies. Standard library only. pytest is a dev
  dependency, not a runtime one.
