#!/usr/bin/env python3

# © 2020 ETH Zurich, Mechanics and Materials Lab
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Displays the differences between two vectors of floats.
These vectors are read from two files.
"""

import math
import argparse
import pathlib
import typing
import sys
import itertools


def float_or_nan(value: str):
    """
    Converts the `value` to float. If this fails then returns NaN.

    >>> float_or_nan("1")
    1.0
    >>> float_or_nan("1.123")
    1.123
    >>> float_or_nan(" 1.123  ")
    1.123
    >>> float_or_nan("ab")
    nan
    >>> float_or_nan(" ")
    nan
    """

    try:
        return float(value)
    except ValueError:
        return math.nan


def diff_iterables(
    sequence_1: typing.Iterable[str],
    sequence_2: typing.Iterable[str],
    rel_tol: float,
    abs_tol: float,
) -> typing.Iterable[str]:
    """
    Reports the differences between the lines of the sequences.

    >>> list(diff_iterables(["1"], ["2"], rel_tol=1, abs_tol=1))
    []
    >>> list(diff_iterables(["1"], ["3"], rel_tol=.7, abs_tol=1))
    []
    >>> list(diff_iterables(["1"], ["3"], rel_tol=.6, abs_tol=1))
    ["Line 0 differs: '1' is not '3'"]
    >>> list(diff_iterables(["1", "-1"], ["3", "3"], rel_tol=.6, abs_tol=1))
    ["Line 0 differs: '1' is not '3'", "Line 1 differs: '-1' is not '3'"]
    >>> list(diff_iterables(["b"], ["a"], rel_tol=1, abs_tol=1))
    []
    >>> list(diff_iterables(["1"], [], rel_tol=1, abs_tol=1))
    ["Line 0 differs: '1' is not ''"]
    >>> list(diff_iterables(["b"], [], rel_tol=1, abs_tol=1))
    []
    """

    for index, strings in enumerate(
        itertools.zip_longest(sequence_1, sequence_2, fillvalue="")
    ):
        values = tuple(float_or_nan(x) for x in strings)

        if not (
            math.isclose(values[0], values[1], rel_tol=rel_tol, abs_tol=abs_tol)
            or (math.isnan(values[0]) and math.isnan(values[1]))
        ):
            yield f"Line {index} differs: '{strings[0].strip()}' is not '{strings[1].strip()}'"


def main():
    """
    Opens the files specified by command line arguments and calls `diff_iterables`.
    Returns the number of differences that were printed.
    """

    parser = argparse.ArgumentParser(
        description="Displays the differences between two vectors of floats. "
        "These vectors are read from two files."
    )
    parser.add_argument(
        "filename_1",
        metavar="file_1",
        type=pathlib.Path,
        help="the file name of the first file",
    )
    parser.add_argument(
        "filename_2",
        metavar="file_2",
        type=pathlib.Path,
        help="the file name of the second file",
    )
    parser.add_argument(
        "--rel_tol", type=float, default=1e-10, help="relative tolerance"
    )
    parser.add_argument(
        "--abs_tol", type=float, default=1e-10, help="absolute tolerance"
    )

    arguments = parser.parse_args()

    with open(str(arguments.filename_1), "r", encoding="utf-8") as file_1:
        with open(str(arguments.filename_2), "r", encoding="utf-8") as file_2:
            return sum(
                1
                for _ in map(
                    print,
                    diff_iterables(
                        file_1,
                        file_2,
                        rel_tol=arguments.rel_tol,
                        abs_tol=arguments.abs_tol,
                    ),
                )
            )


if __name__ == "__main__":
    sys.exit(main())
