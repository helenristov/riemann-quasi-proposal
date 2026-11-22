"""Finite arithmetic experiments for From primes to zeros."""

import argparse
import json

from .arithmetic import at_integer, split_powers
from .poisson import moving_window


def positive_int(arg):
    n = int(arg)
    if n < 1:
        raise argparse.ArgumentTypeError("expected a positive integer")
    return n


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    split = sub.add_parser("split", help="separate first powers, squares and higher exponents")
    split.add_argument("cutoffs", type=positive_int, nargs="+", metavar="X")
    layers = sub.add_parser("layers", help="show the Möbius contributions at each integer")
    layers.add_argument("integers", type=positive_int, nargs="+", metavar="N")
    poisson = sub.add_parser("poisson", help="move a fixed-shape window along the log axis")
    poisson.add_argument("centers", type=positive_int, nargs="+", metavar="CENTER")
    poisson.add_argument("--relative-width", type=float, default=0.125, metavar="R")
    poisson.add_argument("--modes", type=positive_int, default=64, metavar="K")
    args = parser.parse_args()
    try:
        if args.command == "split":
            results = [split_powers(n) for n in args.cutoffs]
        elif args.command == "layers":
            results = [at_integer(n) for n in args.integers]
        else:
            results = [moving_window(c, args.relative_width, args.modes) for c in args.centers]
    except ValueError as exc:
        parser.error(str(exc))
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
