"""Unit tests for `mytools intersect`.

These run without bedtools. Where a test pins something bedtools does that nobody
would derive from first principles, the docstring says so and names the fixture line
that proves it -- the golden suite is what measured it, this is what keeps it.
"""

import argparse
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mytools.bed import Interval, UsageError, overlaps  # noqa: E402
from mytools.intersect import (  # noqa: E402
    BIndex,
    _bin_rank,
    _inflated,
    _region,
    _same_strand,
    run,
)


def iv(chrom, start, end, *rest):
    fields = (chrom, str(start), str(end)) + tuple(rest)
    return Interval(chrom, start, end, fields)


def bed(tmp_path, name, *lines):
    path = tmp_path / name
    path.write_text("".join(line + "\n" for line in lines))
    return str(path)


def intersect(capsys, a, b, **flags):
    """Run the subcommand as the CLI would and return its stdout lines."""
    args = argparse.Namespace(a=a, b=b, u=False, v=False, wa=False, wb=False, s=False)
    for flag, value in flags.items():
        setattr(args, flag, value)
    run(args)
    out = capsys.readouterr().out
    return out.splitlines()


class TestOverlapPredicate:
    """The predicate itself, on the cases issue #7 names.

    `overlaps` is imported from bed.py rather than reimplemented here -- one predicate,
    one place. These pin the behaviour intersect leans on.
    """

    def test_bookended_do_not_overlap(self):
        # a covers 100..199, b starts at 200. They touch and share no base.
        assert not overlaps(iv("chr1", 100, 200), iv("chr1", 200, 300))

    def test_one_base_overlap_does(self):
        # a covers 100..199, b covers 199..298. Base 199 is shared.
        assert overlaps(iv("chr1", 100, 200), iv("chr1", 199, 299))

    def test_fully_nested(self):
        assert overlaps(iv("chr1", 100, 500), iv("chr1", 200, 300))
        assert overlaps(iv("chr1", 200, 300), iv("chr1", 100, 500))

    def test_identical_coordinates(self):
        assert overlaps(iv("chr1", 700, 800), iv("chr1", 700, 800))

    def test_zero_length_bare_predicate(self):
        # Strictly by the predicate: inside an interval it overlaps, sitting on the
        # interval's own start or end it does not. intersect does not use it bare --
        # see TestZeroLengthMatchesOracle.
        assert overlaps(iv("chr1", 100, 200), iv("chr1", 150, 150))
        assert not overlaps(iv("chr1", 100, 200), iv("chr1", 100, 100))
        assert not overlaps(iv("chr1", 100, 200), iv("chr1", 200, 200))

    def test_position_zero(self):
        assert overlaps(iv("chr1", 0, 100), iv("chr1", 0, 1))
        assert not overlaps(iv("chr1", 0, 100), iv("chr1", 100, 200))

    def test_different_chromosomes_never_overlap(self):
        # chr3 is in b.bed and not in a.bed, exactly to exercise this.
        assert not overlaps(iv("chr1", 100, 200), iv("chr3", 100, 200))


class TestInflated:
    def test_zero_length_widens_one_base_each_way(self):
        assert _inflated(iv("chr1", 100, 100))[1:3] == (99, 101)

    def test_zero_length_at_position_zero_widens_below_zero(self):
        # a12 is `chr2 0 0` and bedtools does report it against b10 (chr2 0-10), so
        # the widening is not clamped at 0 -- clamping would leave 0 < 0 false and
        # lose the hit.
        assert _inflated(iv("chr2", 0, 0))[1:3] == (-1, 1)

    def test_ordinary_interval_untouched(self):
        original = iv("chr1", 100, 200)
        assert _inflated(original) is original

    @pytest.mark.parametrize(
        "a,b,expected",
        [
            # Oracle behaviour (bedtools v2.31.1), not reasoning. Each row is a pair
            # from the fixtures and what `intersect` does with it.
            (iv("chr1", 0, 100), iv("chr1", 100, 100), True),    # a01 x b02
            (iv("chr1", 100, 200), iv("chr1", 100, 100), True),  # a02 x b02
            (iv("chr1", 500, 500), iv("chr1", 500, 500), True),  # a07 x b07
            (iv("chr2", 0, 0), iv("chr2", 0, 10), True),         # a12 x b10
            (iv("chr2", 300, 300), iv("chr2", 200, 300), True),  # a16 x b12
            # And where it stops: one base further out and the widening no longer
            # reaches. a17 (chr2 300-400) is bookended with b12 and is in `-v`.
            (iv("chr1", 0, 99), iv("chr1", 100, 100), False),
            (iv("chr1", 101, 200), iv("chr1", 100, 100), False),
            (iv("chr2", 300, 400), iv("chr2", 200, 300), False),
        ],
    )
    def test_zero_length_overlap_matches_oracle(self, a, b, expected):
        assert overlaps(_inflated(a), _inflated(b)) is expected


class TestRegion:
    def test_plain_overlap_is_the_shared_span(self):
        assert _region(iv("chr1", 100, 200), iv("chr1", 150, 250)) == (150, 200)

    def test_nested_b_is_reported_whole(self):
        assert _region(iv("chr1", 100, 500), iv("chr1", 200, 300)) == (200, 300)

    def test_zero_length_b_inflates_the_region(self):
        # Oracle: a01 (chr1 0-100) against b02 (chr1 100-100) prints `chr1 99 100`,
        # and a02 (chr1 100-200) against the same b02 prints `chr1 100 101`. The
        # region is clipped against b widened to 99..101, not against 100..100.
        assert _region(iv("chr1", 0, 100), iv("chr1", 100, 100)) == (99, 100)
        assert _region(iv("chr1", 100, 200), iv("chr1", 100, 100)) == (100, 101)

    def test_zero_length_a_is_not_inflated(self):
        # Oracle, and the asymmetry is the point: a12 (chr2 0-0) against b10
        # (chr2 0-10) prints `chr2 0 0`, not `chr2 0 1`. -a is reported as it was
        # read; only -b widens.
        assert _region(iv("chr2", 0, 0), iv("chr2", 0, 10)) == (0, 0)
        assert _region(iv("chr2", 300, 300), iv("chr2", 200, 300)) == (300, 300)

    def test_both_zero_length(self):
        # a07 x b07: both `chr1 500 500`, reported as `chr1 500 500`.
        assert _region(iv("chr1", 500, 500), iv("chr1", 500, 500)) == (500, 500)


class TestSameStrand:
    def test_same_known_strand_matches(self):
        assert _same_strand(iv("chr1", 0, 1, "x", "0", "+"),
                            iv("chr1", 0, 1, "y", "0", "+"))

    def test_opposite_strands_do_not(self):
        # a03/a04 are the same span on opposite strands; that is what makes -s matter.
        assert not _same_strand(iv("chr1", 150, 250, "a03", "30", "+"),
                                iv("chr1", 150, 250, "a04", "30", "-"))

    def test_unknown_strand_matches_nothing_not_even_itself(self):
        # Oracle: `intersect -s` over two '.'-stranded features that plainly overlap
        # prints nothing, and so does `-s` over two BED3 files. An unknown strand is
        # not a strand -- which is why bed.same_strand() is not used here.
        dot = iv("chr1", 0, 1, "x", "0", ".")
        assert not _same_strand(dot, dot)
        bed3 = iv("chr1", 0, 1)
        assert bed3.strand is None
        assert not _same_strand(bed3, bed3)


class TestBinRank:
    def test_narrow_intervals_share_the_finest_level(self):
        assert _bin_rank(iv("chr1", 100, 200))[0] == 0
        assert _bin_rank(iv("chr1", 100, 100))[0] == 0

    def test_interval_wider_than_a_bin_moves_up_a_level(self):
        # Oracle: a -b feature too wide for the 16kb level is emitted after every
        # narrower hit, whatever order the file put them in.
        assert _bin_rank(iv("chr1", 0, 100000))[0] > _bin_rank(iv("chr1", 100, 200))[0]

    def test_interval_straddling_a_bin_boundary_moves_up(self):
        assert _bin_rank(iv("chr1", 16383, 16385))[0] == 1


class TestBIndex:
    def test_chromosome_absent_from_b_has_no_hits(self, tmp_path):
        # chr3 is in b.bed but not a.bed; the reverse must not blow up either.
        index = BIndex(bed(tmp_path, "b.bed", "chr3\t100\t200\tb19\t27\t+"))
        assert index.hits(iv("chr1", 100, 200)) == []

    def test_hits_come_back_in_b_file_order_not_position_order(self, tmp_path):
        # Oracle: a02 (chr1 100-200) hits b03 (chr1 180-220) *then* b02
        # (chr1 100-100) -- the order they appear in b.bed, because both land in the
        # same finest-level bin. Sorting by start would put b02 first and be wrong.
        path = bed(
            tmp_path, "b.bed",
            "chr1\t180\t220\tb03\t12\t+",
            "chr1\t0\t50\tb01\t11\t+",
            "chr1\t100\t100\tb02\t0\t-",
        )
        names = [b.name for b in BIndex(path).hits(iv("chr1", 100, 200))]
        assert names == ["b03", "b02"]

    def test_wide_b_feature_sorts_after_narrow_ones(self, tmp_path):
        path = bed(
            tmp_path, "b.bed",
            "chr1\t0\t100000\tBIG\t0\t+",
            "chr1\t500\t600\ts1\t0\t+",
            "chr1\t100\t200\ts2\t0\t+",
        )
        names = [b.name for b in BIndex(path).hits(iv("chr1", 0, 1000))]
        assert names == ["s1", "s2", "BIG"]

    def test_strand_filter(self, tmp_path):
        path = bed(
            tmp_path, "b.bed",
            "chr1\t150\t250\tplus\t0\t+",
            "chr1\t150\t250\tminus\t0\t-",
        )
        index = BIndex(path)
        a = iv("chr1", 100, 200, "a", "0", "+")
        assert [b.name for b in index.hits(a)] == ["plus", "minus"]
        assert [b.name for b in index.hits(a, same_strand_only=True)] == ["plus"]

    def test_long_b_feature_is_still_found_from_far_right(self, tmp_path):
        # The prefix-maximum prune must not cut off a feature that starts early and
        # runs long. Without it this is the hit that goes missing.
        lines = ["chr1\t0\t100000\tlong\t0\t+"]
        lines += [f"chr1\t{i}\t{i + 10}\tshort{i}\t0\t+" for i in range(100, 900, 100)]
        index = BIndex(bed(tmp_path, "b.bed", *lines))
        assert [b.name for b in index.hits(iv("chr1", 90000, 90010))] == ["long"]


class TestModes:
    @pytest.fixture
    def files(self, tmp_path):
        a = bed(
            tmp_path, "a.bed",
            "chr1\t100\t200\ta1\t0\t+",   # overlaps b1 and b2
            "chr1\t600\t700\ta2\t0\t-",   # overlaps nothing
            "chr1\t150\t250\ta3\t0\t-",   # overlaps b1 and b2, opposite strand to b1
        )
        b = bed(
            tmp_path, "b.bed",
            "chr1\t180\t220\tb1\t0\t+",
            "chr1\t120\t140\tb2\t0\t-",
        )
        return a, b

    def test_default_prints_the_intersected_region(self, files, capsys):
        a, b = files
        assert intersect(capsys, a, b) == [
            "chr1\t180\t200\ta1\t0\t+",
            "chr1\t120\t140\ta1\t0\t+",
            "chr1\t180\t220\ta3\t0\t-",
        ]

    def test_wa_prints_the_original_a_once_per_hit(self, files, capsys):
        a, b = files
        assert intersect(capsys, a, b, wa=True) == [
            "chr1\t100\t200\ta1\t0\t+",
            "chr1\t100\t200\ta1\t0\t+",
            "chr1\t150\t250\ta3\t0\t-",
        ]

    def test_wb_appends_the_b_feature(self, files, capsys):
        a, b = files
        assert intersect(capsys, a, b, wb=True) == [
            "chr1\t180\t200\ta1\t0\t+\tchr1\t180\t220\tb1\t0\t+",
            "chr1\t120\t140\ta1\t0\t+\tchr1\t120\t140\tb2\t0\t-",
            "chr1\t180\t220\ta3\t0\t-\tchr1\t180\t220\tb1\t0\t+",
        ]

    def test_u_prints_each_a_at_most_once(self, files, capsys):
        a, b = files
        assert intersect(capsys, a, b, u=True) == [
            "chr1\t100\t200\ta1\t0\t+",
            "chr1\t150\t250\ta3\t0\t-",
        ]

    def test_v_prints_the_a_features_with_no_hit(self, files, capsys):
        a, b = files
        assert intersect(capsys, a, b, v=True) == ["chr1\t600\t700\ta2\t0\t-"]

    def test_s_restricts_to_the_same_strand(self, files, capsys):
        a, b = files
        assert intersect(capsys, a, b, s=True) == ["chr1\t180\t200\ta1\t0\t+"]

    def test_input_order_of_a_is_preserved(self, tmp_path, capsys):
        # a.bed is deliberately unsorted and bedtools does not reorder -a.
        a = bed(
            tmp_path, "a.bed",
            "chr2\t100\t200\tsecond\t0\t+",
            "chr1\t100\t200\tfirst\t0\t+",
        )
        b = bed(tmp_path, "b.bed", "chr1\t100\t200\tb\t0\t+", "chr2\t100\t200\tb\t0\t+")
        assert [row.split("\t")[3] for row in intersect(capsys, a, b, u=True)] == [
            "second",
            "first",
        ]

    def test_zero_length_b_at_position_zero(self, tmp_path, capsys):
        """The one place we knowingly differ from the oracle. SPEC.md section 8.

        Real bedtools aborts here -- it bins on `end - 1`, which underflows to -1 --
        so this is our documented deviation, not oracle behaviour, and it is the
        reason this is a unit test and not a golden case. The answer is bedtools'
        own rule carried one base further than bedtools itself can go: a zero-length
        -b feature widens to -1..1, so `chr1 0 10` comes back clipped to `chr1 0 1`,
        exactly as `chr1 5 5` in -b gives `chr1 4 6` where bedtools does run.
        """
        a = bed(tmp_path, "a.bed", "chr1\t0\t10\tq\t0\t+")
        b = bed(tmp_path, "b.bed", "chr1\t0\t0\tz\t0\t+")
        assert intersect(capsys, a, b) == ["chr1\t0\t1\tq\t0\t+"]

        oracle_agrees = bed(tmp_path, "b5.bed", "chr1\t5\t5\tz\t0\t+")
        assert intersect(capsys, a, oracle_agrees) == ["chr1\t4\t6\tq\t0\t+"]

    def test_empty_result_prints_nothing(self, tmp_path, capsys):
        a = bed(tmp_path, "a.bed", "chr1\t100\t200\ta\t0\t+")
        b = bed(tmp_path, "b.bed", "chr9\t100\t200\tb\t0\t+")
        assert intersect(capsys, a, b) == []


class TestUsage:
    @pytest.mark.parametrize(
        "flags",
        [("u", "v"), ("u", "wa"), ("v", "wb"), ("wa", "wb"), ("u", "v", "wa", "wb")],
    )
    def test_mutually_exclusive_modes_are_rejected(self, tmp_path, capsys, flags):
        a = bed(tmp_path, "a.bed", "chr1\t100\t200\ta\t0\t+")
        b = bed(tmp_path, "b.bed", "chr1\t100\t200\tb\t0\t+")
        with pytest.raises(UsageError, match="mutually exclusive"):
            intersect(capsys, a, b, **{flag: True for flag in flags})

    def test_both_inputs_on_stdin_is_rejected(self, capsys):
        with pytest.raises(UsageError, match="stdin"):
            intersect(capsys, "-", "-")
