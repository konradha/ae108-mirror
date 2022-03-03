#!/usr/bin/env python3

# © 2021 ETH Zurich, Mechanics and Materials Lab
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
Runs the tests defined in definition.json files.
"""

import dataclasses
import enum
import itertools
import json
import math
import pathlib
import shutil
import subprocess
import tempfile
import typing
import unittest
import re


ROOT_DIRECTORY = pathlib.Path(__file__).resolve().parent.parent


class ComparisonType(enum.Enum):
    """
    Types of file comparisons.
    """

    NONE = 0
    TEXT = 1
    NUMERIC = 2


@dataclasses.dataclass(frozen=True)
class BinaryOutputDefinition:
    """
    Contains the parameters necessary to compare binary output.
    """

    filename: pathlib.Path
    xdmf_generator_flags: typing.List[str]


@dataclasses.dataclass(frozen=True)
class TestCaseDefinition:
    """
    Contains the parameters necessary to execute a test.
    """

    executable: pathlib.Path
    args: typing.List[str]
    references: pathlib.Path
    compare_stdout: ComparisonType
    ae108_output: BinaryOutputDefinition
    mpi_processes: int


def as_test_case_definitions(
    path: pathlib.Path, definition: typing.Dict[str, typing.Any]
) -> typing.Generator[TestCaseDefinition, None, None]:
    r"""
    Generates test case definitions from `definition` for a test at `path`.

    >>> empty_path = pathlib.Path()
    >>> cwd_path = pathlib.Path.cwd()

    >>> next(
    ...     as_test_case_definitions(
    ...         empty_path,
    ...         {"executable": ["a", "b"]}
    ...     )
    ... ).executable.relative_to(cwd_path)
    PosixPath('a/b')

    >>> next(
    ...     as_test_case_definitions(
    ...         empty_path,
    ...         {"executable": []}
    ...     )
    ... ).args
    []
    >>> next(
    ...     as_test_case_definitions(
    ...         empty_path,
    ...         {"executable": [], "args": ["b"]}
    ...     )
    ... ).args
    ['b']

    >>> next(
    ...     as_test_case_definitions(
    ...         pathlib.Path("a"),
    ...         {"executable": []}
    ...     )
    ... ).references
    PosixPath('a/references')

    >>> next(
    ...     as_test_case_definitions(
    ...         empty_path,
    ...         {"executable": []}
    ...     )
    ... ).compare_stdout
    <ComparisonType.NONE: 0>
    >>> next(
    ...     as_test_case_definitions(
    ...         empty_path,
    ...         {"executable": [], "compare_stdout": "text"}
    ...     )
    ... ).compare_stdout
    <ComparisonType.TEXT: 1>

    >>> next(
    ...     as_test_case_definitions(
    ...         empty_path,
    ...         {"executable": []}
    ...     )
    ... ).mpi_processes
    1
    >>> generator = as_test_case_definitions(
    ...                 empty_path,
    ...                 {"executable": [], "mpi_processes": [1, 2]}
    ...             )
    >>> list(definitions.mpi_processes for definitions in generator)
    [1, 2]
    """
    string_to_comparison_type = {
        "none": ComparisonType.NONE,
        "text": ComparisonType.TEXT,
        "numeric": ComparisonType.NUMERIC,
    }

    for mpi_processes in definition.get("mpi_processes", [1]):
        yield TestCaseDefinition(
            executable=pathlib.Path.cwd() / pathlib.Path(*definition["executable"]),
            args=definition.get("args", []),
            references=path / "references",
            compare_stdout=string_to_comparison_type[
                definition.get("compare_stdout", "none")
            ],
            ae108_output=[
                BinaryOutputDefinition(
                    filename=pathlib.Path(output["filename"]),
                    xdmf_generator_flags=output.get("xdmf_generator_flags", []),
                )
                for output in definition.get("ae108_output", [])
            ],
            mpi_processes=mpi_processes,
        )


def run_executable_with_mpirun(
    executable: pathlib.Path,
    mpi_processes: int,
    args: typing.List[str],
    working_directory: pathlib.Path,
) -> subprocess.CompletedProcess:
    """
    Runs the executable at `executable` with the provided `args`
    from `working_directory` with `mpi_processes` processes.

    >>> empty_path = pathlib.Path()
    >>> run_executable_with_mpirun(
    ...     empty_path,
    ...     mpi_processes = 2,
    ...     args=["-v"],
    ...     working_directory=empty_path
    ... ) # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    subprocess.CalledProcessError: Command '['mpirun', '-n', '2', '.', '-v']' ...
    """
    return run_process(
        args=["mpirun", "-n", str(mpi_processes), str(executable)] + args,
        working_directory=working_directory,
    )


def run_process(
    args: typing.List[str],
    working_directory: pathlib.Path = pathlib.Path.cwd(),
) -> subprocess.CompletedProcess:
    """
    Runs a process with the provided `args` from `working_directory`.

    >>> empty_path = pathlib.Path()
    >>> run_process(
    ...     args=["mpirun", "."],
    ...     working_directory=empty_path
    ... ) # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    subprocess.CalledProcessError: Command '['mpirun', '.']' ...
    """
    return subprocess.run(
        args=args,
        cwd=working_directory,
        capture_output=True,
        check=True,
        text=True,
    )


def diff_text_files(
    case: unittest.TestCase,
    value: pathlib.Path,
    reference: pathlib.Path,
    comparison: ComparisonType,
):
    """
    Compares the files at `value`, `reference` as specified by `comparison.
    Results are reported to `case`.
    Nonexisting references are automatically created.

    >>> path = pathlib.Path(__file__)
    >>> diff_text_files(unittest.TestCase(), path, path, ComparisonType.TEXT)
    """
    if not reference.exists():
        shutil.copy(value, reference)

    comparison_to_function = {
        ComparisonType.TEXT: diff_text_string,
        ComparisonType.NUMERIC: diff_numeric_string,
    }

    with open(value, "r") as value_file:
        with open(reference, "r") as reference_file:
            comparison_to_function.get(comparison, lambda _0, _1, _2: None)(
                value_file.read(), reference_file.read(), case
            )


def diff_vtu_files(
    value: pathlib.Path,
    reference: pathlib.Path,
):
    """
    Compares the files at `value`, `reference`.
    """
    if not reference.exists():
        shutil.copy(value, reference)

    run_process(
        [str(ROOT_DIRECTORY / "tests" / "vtu_diff.py"), str(value), str(reference)]
    )


def diff_text_string(
    value: str, reference: str, case: unittest.TestCase = unittest.TestCase()
):
    """
    Checks that the lines in the strings `value` and `reference` are equal.

    >>> diff_text_string("a", "a")
    >>> diff_text_string("a", "b") # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    AssertionError: 'a' != 'b'
    ...
    """
    case.assertEqual(value, reference)


def float_or_nan(value: str) -> float:
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


def extract_numbers(text: typing.Iterable[str]) -> typing.Iterable[float]:
    """
    Extracts the numbers (separated by ',', '(', ')', ']', '[', '=', or whitespace) from text
    and returns an iterator.

    >>> list(extract_numbers(["1.0, 2.0"]))
    [1.0, 2.0]
    >>> list(extract_numbers(["x=1.0"]))
    [1.0]
    >>> list(extract_numbers(["[1.0]"]))
    [1.0]
    >>> list(extract_numbers(["1.e1"]))
    [10.0]
    >>> list(extract_numbers(["NaN"]))
    []
    >>> list(extract_numbers(["1.0 (2.0)"]))
    [1.0, 2.0]
    >>> list(extract_numbers(["1.0 a 2.0"]))
    [1.0, 2.0]
    >>> list(extract_numbers(["1.0 a 2.0", "b 3.0"]))
    [1.0, 2.0, 3.0]
    """
    return filter(
        lambda x: not math.isnan(x),
        map(
            float_or_nan,
            itertools.chain.from_iterable(
                map(lambda x: re.split(r"[,\(\)\[\]=\s]", x), text)
            ),
        ),
    )


def diff_numeric_string(
    value: str, reference: str, case: unittest.TestCase = unittest.TestCase()
):
    """
    Checks that the lines in the strings `value` and `reference` are almost equal
    when interpreted as floats. Non-float lines are interpreted as NaNs.

    >>> diff_numeric_string("a", "a")
    >>> diff_numeric_string("1", "1")
    >>> diff_numeric_string("1 2", "1 2")
    >>> diff_numeric_string("1", "2")
    Traceback (most recent call last):
    ...
    AssertionError: 1.0 != 2.0 within 7 places (1.0 difference)
    >>> diff_numeric_string("a 1", "b 2")
    Traceback (most recent call last):
    ...
    AssertionError: 1.0 != 2.0 within 7 places (1.0 difference)
    >>> diff_numeric_string("a", "1")
    Traceback (most recent call last):
    ...
    AssertionError: nan != 1.0 within 7 places (nan difference)
    >>> diff_numeric_string("1", "a")
    Traceback (most recent call last):
    ...
    AssertionError: 1.0 != nan within 7 places (nan difference)
    """
    for float_value, float_reference in itertools.zip_longest(
        extract_numbers(iter(value.splitlines())),
        extract_numbers(iter(reference.splitlines())),
        fillvalue=math.nan,
    ):
        case.assertAlmostEqual(float_value, float_reference)


def load_tests(
    loader: unittest.TestLoader,
    standard_tests: unittest.TestSuite,
    _: typing.Any,
) -> unittest.TestSuite:
    """
    Uses the provided definitions in tests/ to create a TestSuite of tests.
    """
    paths = (ROOT_DIRECTORY / "tests").glob("*/*/definition.json")

    for path in paths:
        group_name, test_name = path.parent.parts[-2:]
        to_method_name = (
            lambda processes, name=test_name: f"test_{name}_with_{processes}_mpi_processes"
        )

        with open(path, "r") as file:
            test_case_definitions = as_test_case_definitions(
                path.parent, json.load(file)
            )

            testcase = type(
                group_name,
                (unittest.TestCase,),
                {
                    to_method_name(
                        definition.mpi_processes
                    ): lambda case, definition=definition: run_testcase(
                        definition, case
                    )
                    for definition in test_case_definitions
                },
            )
        standard_tests.addTests(loader.loadTestsFromTestCase(testcase))

    return standard_tests


def run_testcase(
    definition: TestCaseDefinition, case: unittest.TestCase = unittest.TestCase()
):
    """
    Runs the test defined by `definition` and reports the issues to `case`.
    Nonexisting references are automatically created.
    """
    with tempfile.TemporaryDirectory() as directory:
        directory_path = pathlib.Path(directory)

        result = run_executable_with_mpirun(
            executable=definition.executable,
            args=definition.args,
            working_directory=directory_path,
            mpi_processes=definition.mpi_processes,
        )

        with open(directory_path / "stdout.txt", "w") as file:
            file.write(result.stdout)

        diff_text_files(
            case,
            directory_path / "stdout.txt",
            definition.references / "stdout.txt",
            definition.compare_stdout,
        )

        for output in definition.ae108_output:
            with tempfile.TemporaryDirectory() as conversion_directory:
                conversion_path = pathlib.Path(conversion_directory)
                run_process(
                    args=[str(ROOT_DIRECTORY / "tools" / "generate_xdmf.py")]
                    + output.xdmf_generator_flags
                    + [str(directory_path / output.filename)],
                    working_directory=conversion_path,
                )
                run_process(
                    [
                        str(ROOT_DIRECTORY / "tools" / "convert_xdmf_to_vtu.py"),
                        str(conversion_path / output.filename.with_suffix(".xdmf")),
                    ],
                    working_directory=conversion_path,
                )
                diff_vtu_files(
                    conversion_path / output.filename.with_suffix(".vtu"),
                    definition.references / output.filename.with_suffix(".vtu"),
                )


def main() -> None:
    """
    Runs the tests defined in "definition.json"s.
    """
    unittest.main()


if __name__ == "__main__":
    main()
