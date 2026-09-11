"""Argument parsing, dispatch and exit codes.

Exit codes (SPEC.md 7): 0 success, 1 bad input data, 2 bad invocation.
"""

import argparse
import sys

from . import __version__, closest, intersect, merge, sort, subtract
from .bed import BedError, UsageError


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="mytools", description=__doc__)
    parser.add_argument("--version", action="version", version=f"mytools {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="<command>")

    p_sort = subparsers.add_parser("sort", help="sort BED features")
    p_sort.add_argument("-i", required=True, metavar="<file|->")
    p_sort.set_defaults(run=sort.run)

    p_merge = subparsers.add_parser("merge", help="merge overlapping features")
    p_merge.add_argument("-i", required=True, metavar="<file|->")
    p_merge.add_argument("-d", type=int, default=0, metavar="N")
    p_merge.add_argument("-s", action="store_true")
    p_merge.set_defaults(run=merge.run)

    p_intersect = subparsers.add_parser("intersect", help="report overlaps")
    p_intersect.add_argument("-a", required=True, metavar="<file|->")
    p_intersect.add_argument("-b", required=True, metavar="<file>")
    p_intersect.add_argument("-u", action="store_true")
    p_intersect.add_argument("-v", action="store_true")
    p_intersect.add_argument("-wa", action="store_true")
    p_intersect.add_argument("-wb", action="store_true")
    p_intersect.add_argument("-s", action="store_true")
    p_intersect.set_defaults(run=intersect.run)

    p_subtract = subparsers.add_parser("subtract", help="remove -b regions from -a")
    p_subtract.add_argument("-a", required=True, metavar="<file|->")
    p_subtract.add_argument("-b", required=True, metavar="<file>")
    p_subtract.add_argument("-A", action="store_true")
    p_subtract.add_argument("-s", action="store_true")
    p_subtract.set_defaults(run=subtract.run)

    p_closest = subparsers.add_parser("closest", help="find the nearest -b feature")
    p_closest.add_argument("-a", required=True, metavar="<file|->")
    p_closest.add_argument("-b", required=True, metavar="<file>")
    p_closest.add_argument("-d", action="store_true")
    p_closest.add_argument("-io", action="store_true")
    p_closest.add_argument("-s", action="store_true")
    p_closest.set_defaults(run=closest.run)

    return parser


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    parser = build_parser()

    if not argv:
        parser.print_usage(sys.stderr)
        return 2

    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_usage(sys.stderr)
        return 2

    try:
        args.run(args)
    except BrokenPipeError:
        # Downstream closed the pipe (`| head`). Not our problem, and not an error.
        return 0
    except BedError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    except UsageError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    except NotImplementedError:
        print(f"Error: {args.command} is not implemented yet", file=sys.stderr)
        return 1

    return 0
