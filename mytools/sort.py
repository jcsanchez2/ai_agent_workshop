"""mytools sort. See SPEC.md sections 5 and 6.

Sorts BED features by chrom (lexicographic), then start, then end, preserving every
input column verbatim.
"""

import sys

from .bed import read_bed, sort_key


def run(args):
    # sort is the one subcommand that cannot stream: the last line of the input can
    # belong at the top of the output. SPEC.md section 6 targets ~10^6 intervals, which
    # fits comfortably.
    intervals = list(read_bed(args.i))

    # Stable, so features with identical (chrom, start, end) keep their input order --
    # a09/a10 in data/a.bed are byte-identical apart from the name, and bedtools emits
    # them in input order.
    intervals.sort(key=sort_key)

    out = sys.stdout
    for interval in intervals:
        out.write(str(interval))
        out.write("\n")
