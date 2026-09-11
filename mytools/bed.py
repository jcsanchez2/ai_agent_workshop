"""BED parsing and the shared interval logic.

Every subcommand imports from here. If you find yourself writing a second parser or a
second overlap predicate, stop -- that is the bug this module exists to prevent.
"""

import gzip
import io
import sys
from typing import NamedTuple, Optional

GZIP_MAGIC = b"\x1f\x8b"

SKIP_PREFIXES = ("#", "track", "browser")


class BedError(Exception):
    """Bad input data. Exits 1."""


class UsageError(Exception):
    """Bad invocation. Exits 2."""


class Interval(NamedTuple):
    chrom: str
    start: int
    end: int
    fields: tuple

    @property
    def strand(self) -> Optional[str]:
        return self.fields[5] if len(self.fields) > 5 else None

    @property
    def name(self) -> Optional[str]:
        return self.fields[3] if len(self.fields) > 3 else None

    def __str__(self) -> str:
        return "\t".join(self.fields)


def overlaps(a: Interval, b: Interval) -> bool:
    """True if a and b share at least one base.

    BED is 0-based half-open, so `chr1 100 200` covers bases 100..199. The `<` are
    strict on both sides: bookended intervals (a.end == b.start) touch but share no
    base, and do not overlap. Loosening either one to `<=` is the off-by-one this
    whole project is about.
    """
    return a.chrom == b.chrom and a.start < b.end and b.start < a.end


def same_strand(a: Interval, b: Interval) -> bool:
    return a.strand == b.strand


def open_input(path: str):
    """Open a BED source as text, transparently gunzipping it.

    Compression is detected by magic bytes rather than by filename, so a gzipped
    stream arriving on stdin works the same as a .gz file argument.
    """
    if path == "-":
        raw = sys.stdin.buffer
    else:
        try:
            raw = open(path, "rb")
        except FileNotFoundError:
            raise UsageError(f"the requested file ({path}) could not be opened")
        except IsADirectoryError:
            raise UsageError(f"{path} is a directory, not a file")

    if raw.peek(2)[:2] == GZIP_MAGIC:
        raw = gzip.GzipFile(fileobj=raw)

    return io.TextIOWrapper(raw, encoding="utf-8", newline="")


def parse_line(line: str, source: str, lineno: int) -> Optional[Interval]:
    """Parse one BED line. Returns None for lines that should be skipped."""
    stripped = line.rstrip("\n")
    if not stripped.strip():
        return None
    if stripped.startswith(SKIP_PREFIXES):
        return None

    fields = tuple(stripped.split("\t"))
    if len(fields) < 3:
        raise BedError(
            f"{source}:{lineno}: expected at least 3 tab-separated fields, got {len(fields)}"
        )

    try:
        start = int(fields[1])
        end = int(fields[2])
    except ValueError:
        raise BedError(
            f"{source}:{lineno}: non-integer start or end ({fields[1]!r}, {fields[2]!r})"
        )

    if start < 0 or end < 0:
        raise BedError(f"{source}:{lineno}: negative coordinate ({start}, {end})")
    if start > end:
        raise BedError(f"{source}:{lineno}: start > end ({start} > {end})")

    return Interval(fields[0], start, end, fields)


def read_bed(path: str):
    """Yield Intervals from a file path or '-' for stdin."""
    source = "stdin" if path == "-" else path
    with open_input(path) as handle:
        for lineno, line in enumerate(handle, start=1):
            interval = parse_line(line, source, lineno)
            if interval is not None:
                yield interval


def sort_key(interval: Interval):
    """bedtools' default order: chrom lexicographic, then start, then end.

    Lexicographic means chr17 sorts before chr7. That looks wrong and is correct.
    """
    return (interval.chrom, interval.start, interval.end)


def check_sorted(intervals, source: str):
    """Yield intervals, raising if they are not in sort_key order.

    merge and closest require pre-sorted input and error otherwise, same as bedtools.
    """
    previous = None
    for interval in intervals:
        if previous is not None and sort_key(interval) < sort_key(previous):
            raise BedError(
                f"{source}: input is not sorted, out of order record: {interval}"
            )
        previous = interval
        yield interval
