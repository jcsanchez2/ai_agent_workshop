"""mytools merge -- combine overlapping and bookended features.

See SPEC.md 4 (interval semantics), 5 (output) and 6 (memory model).

Three behaviours here are bedtools', not anything you would derive from the spec:
zero-length features are widened by 1bp on each side before merging, `-s` silently
drops features that are not on `+` or `-`, and chromosomes may arrive in any order as
long as each one arrives in a single run. All three were measured against bedtools
v2.31.1 on this VM and are commented where they happen. Do not "fix" them.
"""

import sys
from collections import deque
from itertools import chain, groupby
from operator import attrgetter

from .bed import BedError, check_sorted, read_bed


class _Cluster:
    """One run of features being merged into a single output interval."""

    __slots__ = ("chrom", "start", "end", "count", "first", "closed")

    def __init__(self, interval, start, end):
        self.chrom = interval.chrom
        self.start = start
        self.end = end
        self.count = 1
        self.first = interval
        self.closed = False

    def span(self):
        # A cluster nothing merged into prints the feature's original coordinates:
        # bedtools undoes the zero-length widening when the feature stayed alone, so
        # a lone `chr1 500 500` prints as `500 500` while the same feature beside
        # `chr1 500 600` prints as `499 600`. Oracle behaviour.
        if self.count == 1:
            return self.chrom, self.first.start, self.first.end
        return self.chrom, self.start, self.end


def _widened(interval):
    """The coordinates bedtools merges on.

    A zero-length feature (start == end) is widened to (start - 1, end + 1) first.
    That is why `bedtools merge` on sorted a.bed reports `chr1 499 600` -- the
    zero-length a07 at 500 reaches leftwards to 499 -- and why a zero-length feature
    at position 0 can produce a start of -1. Oracle behaviour; nothing clamps it.
    """
    if interval.start == interval.end:
        return interval.start - 1, interval.end + 1
    return interval.start, interval.end


def merged(intervals, distance=0, stranded=False):
    """Yield (chrom, start, end) for each merged cluster, in bedtools' output order.

    A feature joins the open cluster when the gap between them is <= `distance`.
    That gap test, not bed.overlaps(), is what merge means: bookended features share
    no base -- overlaps() is false for them -- yet their gap is 0, so they merge at
    the default -d 0. Negative `distance` runs the other way and demands that much
    overlap: -d -1 joins only features sharing at least 1bp. SPEC.md 4.

    Clusters come out in the order they were opened, which is how bedtools prints
    them under -s: on `chr1 305 400 -` followed by `chr1 306 306 +, chr1 306 310 +`
    it prints `305 400` before `305 310`, not the lower span first.
    """
    pending = deque()   # clusters in creation order; the front one prints next
    open_clusters = {}  # strand key (None when unstranded) -> cluster still open
    chrom = None

    for interval in intervals:
        if interval.chrom != chrom:
            yield from _close_all(pending, open_clusters)
            chrom = interval.chrom

        key = interval.strand if stranded else None
        if stranded and key not in ("+", "-"):
            # Under -s, bedtools drops features whose strand is "." or absent
            # entirely -- they merge with nothing and never appear. Oracle.
            continue

        start, end = _widened(interval)
        cluster = open_clusters.get(key)
        if cluster is not None and start - cluster.end <= distance:
            cluster.end = max(cluster.end, end)
            cluster.count += 1
        else:
            if cluster is not None:
                cluster.closed = True
            cluster = _Cluster(interval, start, end)
            open_clusters[key] = cluster
            pending.append(cluster)

        yield from _drain(pending)

    yield from _close_all(pending, open_clusters)


def _drain(pending):
    """Emit finished clusters from the front, so creation order is preserved."""
    while pending and pending[0].closed:
        yield pending.popleft().span()


def _close_all(pending, open_clusters):
    """Emit everything still buffered -- end of a chromosome, or end of input."""
    open_clusters.clear()
    while pending:
        yield pending.popleft().span()


def _in_bedtools_order(intervals, source):
    """Yield intervals, rejecting input bedtools would call unsorted (exit 1).

    bedtools does not care about chromosome ORDER -- chr10 before chr2 is accepted --
    it cares that each chromosome arrives in one contiguous run, and that within a run
    starts never go backwards. Ends are not checked: `sort -k1,1 -k2,2n` leaves
    equal-start features in arbitrary end order and bedtools merges them without
    complaint. Measured on v2.31.1.

    bed.check_sorted compares (chrom, start, end), so each run is handed copies with
    end flattened to start; the untouched records are what we yield.
    """
    seen = set()
    for chrom, run in groupby(intervals, key=attrgetter("chrom")):
        run = iter(run)
        first = next(run)
        if chrom in seen:
            raise BedError(
                f"{source}: input is not sorted, out of order record: {first}"
            )
        seen.add(chrom)
        yield from _starts_ascending(chain([first], run), source)


def _starts_ascending(records, source):
    held = deque()

    def flattened():
        for record in records:
            held.append(record)
            yield record._replace(end=record.start)

    for _ in check_sorted(flattened(), source):
        yield held.popleft()


def run(args):
    source = "stdin" if args.i == "-" else args.i
    intervals = _in_bedtools_order(read_bed(args.i), source)

    write = sys.stdout.write
    for chrom, start, end in merged(intervals, args.d, args.s):
        write(f"{chrom}\t{start}\t{end}\n")
