"""Unit tests for `mytools merge`. No bedtools needed; runs in milliseconds.

Cases tagged ORACLE record what bedtools v2.31.1 actually printed on the same input.
They are measurements, not derivations -- if one looks wrong, it is bedtools that is
being strange, and bedtools is the definition of right here. Do not "fix" them.
"""

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from mytools.bed import overlaps, parse_line  # noqa: E402
from mytools.merge import merged  # noqa: E402

MYTOOLS = ROOT / "bin" / "mytools"


def features(*rows):
    """Build Intervals from `(chrom, start, end, strand)` tuples via the real parser."""
    lines = [f"{c}\t{s}\t{e}\tf{n}\t0\t{strand}" for n, (c, s, e, strand) in enumerate(rows)]
    return [parse_line(line, "test", n) for n, line in enumerate(lines, start=1)]


def run_cli(args, stdin=None):
    return subprocess.run(
        [sys.executable, str(MYTOOLS), *args],
        input=stdin,
        capture_output=True,
        text=True,
    )


# --- the gap test -----------------------------------------------------------------

def test_bookended_features_merge_at_d0():
    """Bookended features share no base yet merge at -d 0. Both facts at once.

    overlaps() is false for them -- SPEC.md 4, strict `<` on both sides -- because
    they share no base. merge does not ask about shared bases, it asks whether the
    gap is <= d, and the gap between them is 0.
    """
    a, b = features(("chr1", 100, 200, "+"), ("chr1", 200, 300, "+"))
    assert not overlaps(a, b)
    assert list(merged([a, b])) == [("chr1", 100, 300)]


def test_bookended_features_do_not_merge_at_negative_d():
    """-d -1 means "require 1bp of overlap", and bookended features have none."""
    pair = features(("chr1", 100, 200, "+"), ("chr1", 200, 300, "+"))
    assert list(merged(pair, distance=-1)) == [("chr1", 100, 200), ("chr1", 200, 300)]


def test_one_base_overlap_merges_at_d_minus_one():
    pair = features(("chr1", 100, 200, "+"), ("chr1", 199, 300, "+"))
    assert overlaps(*pair)
    assert list(merged(pair, distance=-1)) == [("chr1", 100, 300)]


def test_one_base_overlap_does_not_merge_at_d_minus_two():
    """-d -2 demands 2bp; a single shared base is not enough."""
    pair = features(("chr1", 100, 200, "+"), ("chr1", 199, 300, "+"))
    assert list(merged(pair, distance=-2)) == [("chr1", 100, 200), ("chr1", 199, 300)]


def test_gap_of_exactly_n_merges_at_d_n():
    pair = features(("chr1", 100, 200, "+"), ("chr1", 210, 300, "+"))
    assert list(merged(pair, distance=10)) == [("chr1", 100, 300)]


def test_gap_of_n_plus_one_does_not_merge_at_d_n():
    pair = features(("chr1", 100, 200, "+"), ("chr1", 211, 300, "+"))
    assert list(merged(pair, distance=10)) == [("chr1", 100, 200), ("chr1", 211, 300)]


def test_nested_features_collapse_to_the_outer_span():
    rows = features(
        ("chr1", 100, 1000, "+"),
        ("chr1", 200, 300, "+"),
        ("chr1", 400, 500, "-"),
    )
    assert list(merged(rows)) == [("chr1", 100, 1000)]


def test_nested_feature_does_not_shorten_the_span():
    """The cluster end is the maximum end seen, not the last one."""
    rows = features(
        ("chr1", 100, 1000, "+"),
        ("chr1", 200, 300, "+"),
        ("chr1", 900, 950, "+"),
    )
    assert list(merged(rows)) == [("chr1", 100, 1000)]


def test_feature_at_position_zero():
    rows = features(("chr1", 0, 100, "+"), ("chr1", 100, 200, "+"))
    assert list(merged(rows)) == [("chr1", 0, 200)]


def test_features_on_different_chromosomes_never_merge():
    rows = features(("chr1", 100, 200, "+"), ("chr2", 150, 250, "+"))
    assert list(merged(rows)) == [("chr1", 100, 200), ("chr2", 150, 250)]


# --- zero-length features: ORACLE ---------------------------------------------------

def test_zero_length_feature_reaches_one_base_leftwards():
    """ORACLE: `bedtools merge` on sorted a.bed reports `chr1 499 600`.

    a07 is the zero-length `chr1 500 500` and a08 is `chr1 500 600`. bedtools widens
    a zero-length feature to (start - 1, end + 1) before merging, so the merged span
    starts at 499. This is the single most surprising thing merge does.
    """
    rows = features(("chr1", 500, 500, "+"), ("chr1", 500, 600, "-"))
    assert list(merged(rows)) == [("chr1", 499, 600)]


def test_lone_zero_length_feature_keeps_its_own_coordinates():
    """ORACLE: the widening is undone when nothing merged with the feature."""
    rows = features(("chr1", 500, 500, "+"), ("chr1", 600, 700, "+"))
    assert list(merged(rows)) == [("chr1", 500, 500), ("chr1", 600, 700)]


def test_zero_length_feature_at_zero_yields_a_negative_start():
    """ORACLE: `chr1 0 0` beside `chr1 0 100` prints `chr1 -1 100`. Nothing clamps."""
    rows = features(("chr1", 0, 0, "+"), ("chr1", 0, 100, "+"))
    assert list(merged(rows)) == [("chr1", -1, 100)]


def test_widening_lets_a_zero_length_feature_bridge_a_gap():
    """ORACLE: `chr1 400 499` and `chr1 500 500` merge at -d 0, into `400 501`."""
    rows = features(("chr1", 400, 499, "+"), ("chr1", 500, 500, "+"))
    assert list(merged(rows)) == [("chr1", 400, 501)]


def test_widening_never_pulls_an_existing_cluster_start_backwards():
    """ORACLE: the span starts where the first feature of the cluster starts."""
    rows = features(("chr1", 100, 200, "+"), ("chr1", 100, 100, "+"))
    assert list(merged(rows)) == [("chr1", 100, 200)]


# --- strand -------------------------------------------------------------------------

def test_stranded_merges_only_within_a_strand():
    rows = features(
        ("chr1", 100, 200, "+"),
        ("chr1", 150, 250, "-"),
        ("chr1", 180, 300, "+"),
    )
    assert list(merged(rows, stranded=True)) == [("chr1", 100, 300), ("chr1", 150, 250)]


def test_stranded_output_is_in_cluster_creation_order_not_coordinate_order():
    """ORACLE: bedtools prints `305 400` before `305 310` on this input."""
    rows = features(
        ("chr1", 305, 400, "-"),
        ("chr1", 306, 306, "+"),
        ("chr1", 306, 310, "+"),
    )
    assert list(merged(rows, stranded=True)) == [("chr1", 305, 400), ("chr1", 305, 310)]


def test_stranded_drops_features_with_no_usable_strand():
    """ORACLE: under -s, `.` (and a missing strand column) vanishes from the output."""
    rows = features(
        ("chr1", 0, 100, "."),
        ("chr1", 100, 200, "+"),
        ("chr1", 200, 300, "."),
    )
    assert list(merged(rows, stranded=True)) == [("chr1", 100, 200)]


def test_unstranded_ignores_strand_entirely():
    rows = features(
        ("chr1", 0, 100, "."),
        ("chr1", 100, 200, "+"),
        ("chr1", 200, 300, "-"),
    )
    assert list(merged(rows)) == [("chr1", 0, 300)]


def test_empty_input():
    assert list(merged([])) == []


# --- the command line ----------------------------------------------------------------

def test_output_is_bed3_with_input_columns_dropped():
    result = run_cli(["merge", "-i", "-"], stdin="chr1\t100\t200\tname\t60\t+\n")
    assert result.returncode == 0
    assert result.stdout == "chr1\t100\t200\n"


def test_reads_from_stdin():
    text = "chr1\t0\t100\ta\t0\t+\nchr1\t100\t200\tb\t0\t-\n"
    result = run_cli(["merge", "-i", "-"], stdin=text)
    assert result.returncode == 0
    assert result.stdout == "chr1\t0\t200\n"


def test_unsorted_input_exits_1_and_names_the_offending_line():
    text = "chr1\t300\t400\ta\t0\t+\nchr1\t100\t200\tb\t0\t+\n"
    result = run_cli(["merge", "-i", "-"], stdin=text)
    assert result.returncode == 1
    assert result.stdout == ""
    assert "chr1\t100\t200\tb\t0\t+" in result.stderr


def test_revisited_chromosome_exits_1():
    text = "chr1\t0\t10\ta\t0\t+\nchr2\t0\t10\tb\t0\t+\nchr1\t0\t10\tc\t0\t+\n"
    result = run_cli(["merge", "-i", "-"], stdin=text)
    assert result.returncode == 1
    assert "chr1\t0\t10\tc\t0\t+" in result.stderr


def test_chromosomes_may_arrive_in_any_order():
    """ORACLE: bedtools accepts chr10 before chr2 -- only revisits are an error."""
    text = "chr10\t100\t300\ta\t0\t+\nchr2\t100\t200\tb\t0\t+\n"
    result = run_cli(["merge", "-i", "-"], stdin=text)
    assert result.returncode == 0
    assert result.stdout == "chr10\t100\t300\nchr2\t100\t200\n"


def test_equal_starts_with_decreasing_ends_are_not_unsorted():
    """ORACLE: `sort -k1,1 -k2,2n` leaves these ties unordered and bedtools accepts it."""
    text = "chr1\t100\t300\ta\t0\t+\nchr1\t100\t200\tb\t0\t+\n"
    result = run_cli(["merge", "-i", "-"], stdin=text)
    assert result.returncode == 0
    assert result.stdout == "chr1\t100\t300\n"


def test_missing_input_file_is_a_usage_error():
    result = run_cli(["merge", "-i", "no/such/file.bed"])
    assert result.returncode == 2


@pytest.mark.parametrize("flag", ["-d", "-s"])
def test_flags_are_accepted(flag):
    args = ["merge", flag, "1", "-i", "-"] if flag == "-d" else ["merge", flag, "-i", "-"]
    assert run_cli(args, stdin="chr1\t0\t10\ta\t0\t+\n").returncode == 0


def test_negative_distance_on_the_command_line():
    text = "chr1\t100\t200\ta\t0\t+\nchr1\t200\t300\tb\t0\t+\n"
    result = run_cli(["merge", "-d", "-1", "-i", "-"], stdin=text)
    assert result.returncode == 0
    assert result.stdout == "chr1\t100\t200\nchr1\t200\t300\n"
