"""Unit tests for `mytools sort`.

These run without bedtools. Each one pins a single behaviour that the golden case
(tests/run_golden.sh, SPEC.md section 8 case 1) only covers in aggregate: when the
golden diff fails, one of these should name the reason.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mytools.cli import main  # noqa: E402


def run_sort(capsys, path):
    """Run `mytools sort -i <path>` in-process. Returns (exit code, stdout lines)."""
    rc = main(["sort", "-i", str(path)])
    out = capsys.readouterr().out
    return rc, out.splitlines()


def write(tmp_path, text, name="in.bed"):
    path = tmp_path / name
    path.write_text(text)
    return path


def names(lines):
    return [line.split("\t")[3] for line in lines]


class TestOrder:
    def test_lexicographic_chrom_order(self, tmp_path, capsys):
        # chr17 before chr2 before chr7 before chrX. Lexicographic, not numeric:
        # it looks wrong to humans and it is bedtools' default.
        path = write(
            tmp_path,
            "chr7\t0\t1\tg\t0\t+\n"
            "chrX\t0\t1\tx\t0\t+\n"
            "chr17\t0\t1\tq\t0\t+\n"
            "chr2\t0\t1\tb\t0\t+\n",
        )
        rc, lines = run_sort(capsys, path)
        assert rc == 0
        assert names(lines) == ["q", "b", "g", "x"]

    def test_start_orders_within_a_chrom(self, tmp_path, capsys):
        path = write(
            tmp_path,
            "chr1\t300\t400\tthird\t0\t+\n"
            "chr1\t100\t500\tfirst\t0\t+\n"
            "chr1\t200\t250\tsecond\t0\t+\n",
        )
        rc, lines = run_sort(capsys, path)
        assert rc == 0
        assert names(lines) == ["first", "second", "third"]

    def test_end_breaks_ties_on_start(self, tmp_path, capsys):
        path = write(
            tmp_path,
            "chr1\t100\t300\twide\t0\t+\n"
            "chr1\t100\t200\tnarrow\t0\t+\n"
            "chr1\t100\t100\tempty\t0\t+\n",
        )
        rc, lines = run_sort(capsys, path)
        assert rc == 0
        assert names(lines) == ["empty", "narrow", "wide"]

    def test_identical_intervals_keep_input_order(self, tmp_path, capsys):
        # a09/a10 in data/a.bed differ only by name, and bedtools emits them in
        # input order. The sort must be stable.
        path = write(
            tmp_path,
            "chr1\t700\t800\ta09\t60\t+\n"
            "chr1\t700\t800\ta10\t60\t+\n",
        )
        rc, lines = run_sort(capsys, path)
        assert rc == 0
        assert names(lines) == ["a09", "a10"]

    def test_position_zero_sorts_first(self, tmp_path, capsys):
        # Position 0 is a real coordinate, not a missing value.
        path = write(
            tmp_path,
            "chr1\t100\t200\tlater\t0\t+\n"
            "chr1\t0\t100\tat_zero\t0\t+\n"
            "chr1\t0\t0\tzero_length_at_zero\t0\t+\n",
        )
        rc, lines = run_sort(capsys, path)
        assert rc == 0
        assert names(lines) == ["zero_length_at_zero", "at_zero", "later"]


class TestColumnsPreserved:
    def test_bed12_survives_with_all_columns(self, tmp_path, capsys):
        line = "\t".join(
            [
                "chr1", "1000", "2000", "feature", "960", "+", "1000", "2000",
                "255,0,0", "2", "100,200", "0,800",
            ]
        )
        path = write(tmp_path, line + "\n")
        rc, lines = run_sort(capsys, path)
        assert rc == 0
        assert lines == [line]
        assert len(lines[0].split("\t")) == 12

    def test_ragged_column_counts_are_each_preserved(self, tmp_path, capsys):
        # SPEC.md section 3: column count may vary between lines, and trailing
        # columns are carried through untouched rather than parsed.
        path = write(
            tmp_path,
            "chr1\t300\t400\tnamed\t0\t+\n"
            "chr1\t100\t200\n",
        )
        rc, lines = run_sort(capsys, path)
        assert rc == 0
        assert lines == ["chr1\t100\t200", "chr1\t300\t400\tnamed\t0\t+"]

    def test_comments_and_track_lines_are_dropped(self, tmp_path, capsys):
        path = write(
            tmp_path,
            "# a comment\n"
            "track name=stuff\n"
            "chr1\t100\t200\tonly\t0\t+\n"
            "\n",
        )
        rc, lines = run_sort(capsys, path)
        assert rc == 0
        assert lines == ["chr1\t100\t200\tonly\t0\t+"]

    def test_empty_input_prints_nothing_and_succeeds(self, tmp_path, capsys):
        rc, lines = run_sort(capsys, write(tmp_path, ""))
        assert rc == 0
        assert lines == []


class TestExitCodes:
    """SPEC.md section 7: 1 for bad data, 2 for bad invocation."""

    @pytest.mark.parametrize(
        "line",
        [
            "chr1\tfoo\t200\n",   # non-integer coordinate
            "chr1\t500\t400\n",   # start > end
            "chr1\t100\n",        # fewer than 3 fields
        ],
    )
    def test_malformed_data_exits_1(self, tmp_path, capsys, line):
        rc = main(["sort", "-i", str(write(tmp_path, line))])
        assert rc == 1
        assert capsys.readouterr().out == ""

    def test_errors_go_to_stderr_not_stdout(self, tmp_path, capsys):
        main(["sort", "-i", str(write(tmp_path, "chr1\t500\t400\n"))])
        captured = capsys.readouterr()
        assert captured.out == ""
        assert "start > end" in captured.err

    def test_missing_file_exits_2(self, tmp_path, capsys):
        # A deliberate deviation: bedtools exits 1 here. SPEC.md section 7 explains
        # why, and confines it to cases no golden test adjudicates.
        rc = main(["sort", "-i", str(tmp_path / "nope.bed")])
        assert rc == 2
        assert capsys.readouterr().out == ""
