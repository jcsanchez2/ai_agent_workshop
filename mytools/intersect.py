"""mytools intersect -- report overlaps between -a and -b. See SPEC.md 4, 5, 6.

`-b` is loaded into memory grouped by chrom and sorted by start, then binary-searched;
`-a` is streamed and its input order preserved.

Two behaviours here look wrong and are not. Both are measured against the oracle
(bedtools v2.31.1) and are explained where they happen:

  * zero-length intervals are inflated 1bp each way before the overlap test, and the
    inflated `-b` interval is what clips the reported region (`_inflated`, `_region`);
  * when one `-a` feature hits several `-b` features, bedtools emits them in the order
    of its bin index, not by position (`_bin_rank`).
"""

import sys
from bisect import bisect_left, bisect_right
from collections import defaultdict

from .bed import UsageError, overlaps, read_bed

#: The four output modes, in the order the usage message lists them. Mutually exclusive.
MODES = ("u", "v", "wa", "wb")

#: bedtools only counts a strand as matching when it is actually known.
KNOWN_STRANDS = ("+", "-")


def _inflated(interval):
    """A zero-length interval widened 1bp each way; anything else unchanged.

    Oracle behaviour, not reasoning. `overlaps()` says a zero-length interval shares no
    base with anything whose end it sits on, but bedtools disagrees: `b02` is
    `chr1 100 100` and `intersect -a data/a.bed -b data/b.bed` reports it against both
    `a01` (chr1 0-100) and `a02` (chr1 100-200). Widening to `chr1 99 101` reproduces
    every case in the fixtures, on both sides of the comparison -- `a07`/`a12`/`a16` are
    zero-length in `-a` and behave the same way.

    One case the oracle cannot adjudicate: a zero-length interval at position 0 in
    `-b` widens to `-1..1`, and bedtools aborts rather than handling it. We do not
    follow it there -- see SPEC.md section 8, which is where that deviation is recorded.

    The widened form is used for the overlap test and for clipping the reported region
    only. Printing always uses the interval as it was read, which is why `-wb` shows
    `chr1 100 100` for `b02` while the region beside it was clipped against `99 101`.
    """
    if interval.start == interval.end:
        return interval._replace(start=interval.start - 1, end=interval.end + 1)
    return interval


def _same_strand(a, b):
    """bedtools' `-s` rule: both strands known, and equal.

    Deliberately not `bed.same_strand()`, which compares the raw values and so would
    call two `.`-stranded features (or two BED3 features, which have no strand column
    at all) a match. bedtools treats an unknown strand as matching nothing: `intersect
    -s` over two BED3 files, or over two `.`-stranded features that plainly overlap,
    prints nothing. Measured on v2.31.1.
    """
    return a.strand in KNOWN_STRANDS and a.strand == b.strand


_BIN_FIRST_SHIFT = 14  # finest bin is 2**14 == 16kb
_BIN_NEXT_SHIFT = 3    # each level up is 8x wider
_BIN_LEVELS = 7


def _bin_rank(interval):
    """Sort key reproducing the order bedtools emits hits in: (bin level, bin, input).

    bedtools indexes `-b` in a UCSC-style bin hierarchy -- 7 levels, the finest 16kb
    wide, each level 8x the one below -- and answers a query by walking the levels
    finest-first, bins ascending, features within a bin in the order they were read.
    That, not position, is the order hits come out in, and `intersect` prints them in
    the order it gets them.

    `a02` (chr1 100-200) is the case that proves it: its hits are `b03` (chr1 180-220)
    then `b02` (chr1 100-100), which is neither sorted by start nor sorted by name --
    it is the order the two appear in `data/b.bed`, because both land in the same
    finest-level bin. Mixing in a feature too wide for that bin moves it after every
    narrower hit regardless of where the file put it, which is what the level term is
    for.
    """
    # Half-open: the last base is end-1. A zero-length interval has no last base, so
    # it is binned on its start.
    last_base = interval.end - 1 if interval.end > interval.start else interval.start
    lo = interval.start >> _BIN_FIRST_SHIFT
    hi = last_base >> _BIN_FIRST_SHIFT
    for level in range(_BIN_LEVELS):
        if lo == hi:
            return (level, lo)
        lo >>= _BIN_NEXT_SHIFT
        hi >>= _BIN_NEXT_SHIFT
    # Wider than the whole hierarchy. bedtools drops these in the single top bin.
    return (_BIN_LEVELS, 0)


class BIndex:
    """`-b`, held in memory: grouped by chrom, sorted by start, binary-searched.

    SPEC.md 6. `-a` streams past this; nothing here is quadratic in the size of `-b`.
    Coordinates are stored inflated (see `_inflated`) because that is what the overlap
    test uses; the Interval kept alongside is the one that was read, and is what gets
    printed.
    """

    def __init__(self, path):
        grouped = defaultdict(list)
        for order, interval in enumerate(read_bed(path)):
            grouped[interval.chrom].append((_inflated(interval), interval, order))

        self._by_chrom = {}
        for chrom, rows in grouped.items():
            rows.sort(key=lambda row: (row[0].start, row[0].end))
            wide = [row[0] for row in rows]
            starts = [row[0].start for row in rows]
            intervals = [row[1] for row in rows]
            ranks = [_bin_rank(row[1]) + (row[2],) for row in rows]

            # Running maximum of the ends seen so far. It only goes up, so it can be
            # binary-searched for the leftmost feature that can still reach the query
            # -- without it, a query would have to scan every feature that starts
            # before it, which is the quadratic failure SPEC.md 6 is about.
            reach = []
            highest = -1
            for interval in wide:
                highest = max(highest, interval.end)
                reach.append(highest)

            self._by_chrom[chrom] = (starts, wide, reach, intervals, ranks)

    def hits(self, a, same_strand_only=False):
        """Every `-b` feature overlapping `a`, in the order bedtools emits them."""
        entry = self._by_chrom.get(a.chrom)
        if entry is None:
            return []
        starts, wide, reach, intervals, ranks = entry

        query = _inflated(a)
        # Narrow to the features that could possibly overlap before testing any of
        # them: those starting before the query ends, and reaching past where it
        # starts. `overlaps()` still decides -- this only bounds how often it is asked.
        last = bisect_left(starts, query.end)
        first = bisect_right(reach, query.start)

        found = []
        for i in range(first, last):
            if not overlaps(query, wide[i]):
                continue  # inside the window, but stops short of the query
            b = intervals[i]
            if same_strand_only and not _same_strand(a, b):
                continue
            found.append((ranks[i], b))

        found.sort(key=lambda row: row[0])
        return [b for _, b in found]


def _region(a, b):
    """The intersected region of `a` and `b`, as bedtools reports it.

    The `-b` side is clipped against its inflated form and the `-a` side against its
    literal one -- bedtools is asymmetric here. `a12` (chr2 0-0) against `b10`
    (chr2 0-10) reports `chr2 0 0`, so `-a` is not widened; `a01` (chr1 0-100) against
    `b02` (chr1 100-100) reports `chr1 99 100`, so `-b` is.
    """
    wide_b = _inflated(b)
    return max(a.start, wide_b.start), min(a.end, wide_b.end)


def _line(*parts):
    return "\t".join(parts) + "\n"


def _region_fields(a, b):
    """The default output row: the intersected region, carrying `-a`'s extra columns."""
    start, end = _region(a, b)
    return (a.chrom, str(start), str(end)) + a.fields[3:]


def run(args):
    modes = [mode for mode in MODES if getattr(args, mode)]
    if len(modes) > 1:
        raise UsageError(
            "-{} are mutually exclusive; pick one".format(", -".join(modes))
        )
    mode = modes[0] if modes else None

    if args.a == "-" and args.b == "-":
        raise UsageError("only one of -a and -b may be stdin")

    index = BIndex(args.b)
    out = sys.stdout

    for a in read_bed(args.a):
        hits = index.hits(a, same_strand_only=args.s)

        if mode == "v":
            if not hits:
                out.write(_line(*a.fields))
        elif mode == "u":
            if hits:
                out.write(_line(*a.fields))
        elif mode == "wa":
            for _ in hits:
                out.write(_line(*a.fields))
        elif mode == "wb":
            for b in hits:
                out.write(_line(*_region_fields(a, b), *b.fields))
        else:
            for b in hits:
                out.write(_line(*_region_fields(a, b)))
