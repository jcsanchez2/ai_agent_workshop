"""Unit tests for the shared interval logic.

These run without bedtools. Each one pins exactly one behaviour, so a failure names
the bug rather than just reporting disagreement.
"""

import gzip
import io
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mytools.bed import (  # noqa: E402
    BedError,
    Interval,
    check_sorted,
    overlaps,
    parse_line,
    read_bed,
    sort_key,
)


def iv(chrom, start, end, *rest):
    fields = (chrom, str(start), str(end)) + tuple(rest)
    return Interval(chrom, start, end, fields)


class TestOverlaps:
    def test_plain_overlap(self):
        assert overlaps(iv("chr1", 100, 200), iv("chr1", 150, 250))

    def test_one_base_overlap(self):
        # a covers 100..198, b covers 198..299. They share base 198.
        assert overlaps(iv("chr1", 100, 199), iv("chr1", 198, 300))

    def test_bookended_do_not_overlap(self):
        # a covers 100..199, b starts at 200. They touch but share no base.
        assert not overlaps(iv("chr1", 100, 200), iv("chr1", 200, 300))

    def test_nested(self):
        assert overlaps(iv("chr1", 100, 500), iv("chr1", 200, 300))

    def test_identical(self):
        assert overlaps(iv("chr1", 100, 200), iv("chr1", 100, 200))

    def test_different_chromosomes_never_overlap(self):
        assert not overlaps(iv("chr1", 100, 200), iv("chr2", 100, 200))

    def test_position_zero(self):
        assert overlaps(iv("chr1", 0, 100), iv("chr1", 0, 50))

    def test_zero_length_inside_interval_overlaps(self):
        # Intuition says a zero-length interval has no bases and so overlaps nothing.
        # The predicate says otherwise (100 < 150 and 150 < 200), and bedtools agrees:
        # `intersect -u` reports the feature. Verified against the oracle, because
        # guessing here is how you enshrine your own misunderstanding.
        #
        # What bedtools does with the *reported region* is separate and stranger -- it
        # prints chr1 149 151, inflating the zero-length interval one base each way.
        # That belongs to intersect (#7) and its golden tests, not to the predicate.
        assert overlaps(iv("chr1", 100, 200), iv("chr1", 150, 150))

    def test_zero_length_at_interval_start_does_not_overlap(self):
        # b.start == a.start == 100, b.end == 100, so a.start < b.end is false.
        assert not overlaps(iv("chr1", 100, 200), iv("chr1", 100, 100))

    def test_disjoint(self):
        assert not overlaps(iv("chr1", 100, 200), iv("chr1", 300, 400))


class TestSortKey:
    def test_lexicographic_chrom_order(self):
        # chr17 before chr7 -- lexicographic, not numeric. Looks wrong, is correct.
        assert sort_key(iv("chr17", 0, 1)) < sort_key(iv("chr7", 0, 1))

    def test_start_before_end(self):
        assert sort_key(iv("chr1", 100, 500)) < sort_key(iv("chr1", 200, 300))

    def test_end_breaks_start_ties(self):
        assert sort_key(iv("chr1", 100, 200)) < sort_key(iv("chr1", 100, 300))


class TestParseLine:
    def test_bed3(self):
        assert parse_line("chr1\t100\t200", "x", 1) == iv("chr1", 100, 200)

    def test_bed6_keeps_strand(self):
        assert parse_line("chr1\t100\t200\tn\t0\t-", "x", 1).strand == "-"

    def test_bed12_preserves_all_columns(self):
        line = "\t".join(["chr1", "100", "200"] + [str(i) for i in range(9)])
        assert len(parse_line(line, "x", 1).fields) == 12

    def test_zero_length_is_legal(self):
        assert parse_line("chr1\t100\t100", "x", 1) == iv("chr1", 100, 100)

    def test_position_zero_is_legal(self):
        assert parse_line("chr1\t0\t0", "x", 1) == iv("chr1", 0, 0)

    @pytest.mark.parametrize("line", ["", "   ", "# comment", "track name=x", "browser x"])
    def test_skipped_lines(self, line):
        assert parse_line(line, "x", 1) is None

    def test_start_greater_than_end_is_an_error(self):
        with pytest.raises(BedError, match="start > end"):
            parse_line("chr1\t500\t400", "a.bed", 14)

    def test_error_names_file_and_line(self):
        with pytest.raises(BedError, match=r"a\.bed:14"):
            parse_line("chr1\t500\t400", "a.bed", 14)

    def test_non_integer_coordinate(self):
        with pytest.raises(BedError, match="non-integer"):
            parse_line("chr1\tfoo\t200", "x", 1)

    def test_too_few_fields(self):
        with pytest.raises(BedError, match="at least 3"):
            parse_line("chr1\t100", "x", 1)


class TestCheckSorted:
    def test_sorted_passes(self):
        rows = [iv("chr1", 0, 10), iv("chr1", 10, 20), iv("chr2", 0, 5)]
        assert list(check_sorted(rows, "x")) == rows

    def test_unsorted_raises(self):
        rows = [iv("chr1", 100, 200), iv("chr1", 0, 50)]
        with pytest.raises(BedError, match="not sorted"):
            list(check_sorted(rows, "x"))


class TestReadBed:
    def test_reads_plain_file(self, tmp_path):
        p = tmp_path / "x.bed"
        p.write_text("chr1\t100\t200\nchr1\t300\t400\n")
        assert len(list(read_bed(str(p)))) == 2

    def test_gzip_detected_by_magic_bytes_not_extension(self, tmp_path):
        # Deliberately not named .gz -- detection must not rely on the filename.
        p = tmp_path / "x.bed"
        p.write_bytes(gzip.compress(b"chr1\t100\t200\n"))
        assert list(read_bed(str(p))) == [iv("chr1", 100, 200)]

    def test_stdin(self, monkeypatch):
        buf = io.BytesIO(b"chr1\t100\t200\n")
        monkeypatch.setattr(sys, "stdin", io.TextIOWrapper(io.BufferedReader(buf)))
        assert list(read_bed("-")) == [iv("chr1", 100, 200)]
