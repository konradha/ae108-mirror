#!/usr/bin/env python3

# © 2022 ETH Zurich, Mechanics and Materials Lab
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
Compares the contents of two VTU files.
"""

import unittest
import typing
import pathlib
import argparse
import sys
from functools import partial
import vtk
import scipy.spatial

RESULT_GRID = vtk.vtkUnstructuredGrid()
REFERENCE_GRID = vtk.vtkUnstructuredGrid()


def point_coordinates(
    grid: vtk.vtkUnstructuredGrid, point_id: int
) -> typing.List[float]:
    """
    Returns the coordinates of the cell with index `point_id`.
    """
    return list(grid.GetPoint(point_id))


def all_point_coordinates(
    grid: vtk.vtkUnstructuredGrid,
) -> typing.List[typing.List[float]]:
    """
    Returns a list of point coordinates.
    """
    return list(map(partial(point_coordinates, grid), range(grid.GetNumberOfPoints())))


def point_index_permutation() -> typing.List[int]:
    """
    Infers the result index of a point based on its coordinates.
    """
    result = all_point_coordinates(RESULT_GRID)
    reference = all_point_coordinates(REFERENCE_GRID)

    tree = scipy.spatial.KDTree(result)
    return list(tree.query(reference)[1])


def cell_connectivity(grid: vtk.vtkUnstructuredGrid, cell_id: int) -> typing.List[int]:
    """
    Returns the connectivity for the cell with index `cell_id`.
    """
    ids = vtk.vtkIdList()
    grid.GetCells().GetCellAtId(cell_id, ids)
    return list(map(ids.GetId, range(ids.GetNumberOfIds())))


def all_cell_connectivities(
    grid: vtk.vtkUnstructuredGrid,
) -> typing.List[typing.List[int]]:
    """
    Returns the connectivity for all cells in `grid`.
    """
    return list(map(partial(cell_connectivity, grid), range(grid.GetNumberOfCells())))


def cell_index_permutation() -> typing.List[int]:
    """
    Infers the result index of a cell based on its connectivity.
    """
    point_permutation = point_index_permutation()
    reference = map(
        lambda c: tuple(map(point_permutation.__getitem__, c)),
        all_cell_connectivities(REFERENCE_GRID),
    )
    result = {tuple(t): i for i, t in enumerate(all_cell_connectivities(RESULT_GRID))}

    return list(map(result.__getitem__, reference))


DataType = typing.TypeVar("DataType")


def permute_point_data(data: typing.List[DataType]) -> typing.Iterable[DataType]:
    """
    Permutes the items of `data` into the reference order inferred from the coordinates.
    """
    permutation = point_index_permutation()
    for i in permutation:
        yield data[i]


def permute_cell_data(data: typing.List[DataType]) -> typing.Iterable[DataType]:
    """
    Permutes the items of `data` into the reference order inferred from the connectivities.
    """
    permutation = cell_index_permutation()
    for i in permutation:
        yield data[i]


def all_point_data(
    grid: vtk.vtkUnstructuredGrid,
) -> typing.List[typing.Dict[str, typing.Any]]:
    """
    Returns a list of real-valued point data.
    Ignores arrays that end with ".imag" and drops a ".real" suffix in the name.
    """
    data = grid.GetPointData()
    return [
        {
            "name": data.GetArrayName(i)
            if not data.GetArrayName(i).endswith(".real")
            else data.GetArrayName(i)[:-5],
            "array": [
                list(data.GetArray(i).GetTuple(j))
                for j in range(data.GetArray(i).GetNumberOfTuples())
            ],
        }
        for i in range(data.GetNumberOfArrays())
        if not data.GetArrayName(i).endswith(".imag")
    ]


def all_cell_data(
    grid: vtk.vtkUnstructuredGrid,
) -> typing.List[typing.Dict[str, typing.Any]]:
    """
    Returns a list of real-valued cell data.
    Ignores arrays that end with ".imag" and drops a ".real" suffix in the name.
    """
    data = grid.GetCellData()
    return [
        {
            "name": data.GetArrayName(i)
            if not data.GetArrayName(i).endswith(".real")
            else data.GetArrayName(i)[:-5],
            "array": [
                list(data.GetArray(i).GetTuple(j))
                for j in range(data.GetArray(i).GetNumberOfTuples())
            ],
        }
        for i in range(data.GetNumberOfArrays())
        if not data.GetArrayName(i).endswith(".imag")
    ]


def all_cell_types(
    grid: vtk.vtkUnstructuredGrid,
) -> typing.List[int]:
    """
    Returns a list of all cell types.
    """
    return list(map(grid.GetCellType, range(grid.GetNumberOfCells())))


class GridComparison(unittest.TestCase):
    """
    Compares two grids.
    """

    def test_same_number_of_points(self):
        """
        The grids have the same number of points.
        """
        self.assertEqual(
            RESULT_GRID.GetNumberOfPoints(), REFERENCE_GRID.GetNumberOfPoints()
        )

    def test_same_number_of_cells(self):
        """
        The grids have the same number of cells.
        """
        self.assertEqual(
            RESULT_GRID.GetNumberOfCells(), REFERENCE_GRID.GetNumberOfCells()
        )

    def test_same_point_coordinates(self):
        """
        The points of the grids have the same coordinates.
        """
        for result, reference in zip(
            permute_point_data(all_point_coordinates(RESULT_GRID)),
            all_point_coordinates(REFERENCE_GRID),
        ):
            for result_x, reference_x in zip(result, reference):
                self.assertAlmostEqual(result_x, reference_x)

    def test_same_cell_connectivities(self):
        """
        The cells of the grids have the same vertex indices (up to permutation).
        """
        permutation = point_index_permutation()
        for result, reference in zip(
            permute_cell_data(all_cell_connectivities(RESULT_GRID)),
            all_cell_connectivities(REFERENCE_GRID),
        ):
            for result_x, reference_x in zip(
                result, map(permutation.__getitem__, reference)
            ):
                self.assertAlmostEqual(result_x, reference_x)

    def test_same_cell_types(self):
        """
        The cells of the grid have the same types.
        """
        for result, reference in zip(
            permute_cell_data(all_cell_types(RESULT_GRID)),
            all_cell_types(REFERENCE_GRID),
        ):
            self.assertEqual(result, reference)

    def test_same_number_of_real_point_data(self):
        """
        The points of the grid have the same number of data.
        """
        self.assertEqual(
            len(all_point_data(RESULT_GRID)), len(all_point_data(REFERENCE_GRID))
        )

    def test_same_real_point_data_arrays(self):
        """
        The points of the grid have the same point data at the same coordinates.
        """
        for result, reference in zip(
            sorted(all_point_data(RESULT_GRID), key=lambda x: x["name"]),
            sorted(all_point_data(REFERENCE_GRID), key=lambda x: x["name"]),
        ):
            self.assertEqual(result["name"], reference["name"])
            for result_item, reference_item in zip(
                permute_point_data(result["array"]), reference["array"]
            ):
                for result_x, reference_x in zip(result_item, reference_item):
                    self.assertAlmostEqual(result_x, reference_x)

    def test_same_number_of_real_cell_data(self):
        """
        The cells of the grid have the same number of data.
        """
        self.assertEqual(
            len(all_cell_data(RESULT_GRID)), len(all_cell_data(REFERENCE_GRID))
        )

    def test_same_real_cell_data_arrays(self):
        """
        The cells of the grid have the same cell data.
        """

        for result, reference in zip(
            sorted(all_cell_data(RESULT_GRID), key=lambda x: x["name"]),
            sorted(all_cell_data(REFERENCE_GRID), key=lambda x: x["name"]),
        ):
            self.assertEqual(result["name"], reference["name"])
            for result_item, reference_item in zip(
                permute_cell_data(result["array"]), reference["array"]
            ):
                for result_x, reference_x in zip(result_item, reference_item):
                    self.assertAlmostEqual(result_x, reference_x)


def read_vtu(path: pathlib.Path) -> vtk.vtkUnstructuredGrid:
    """
    Reads the VTU file at path and returns the result.
    """
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    return reader.GetOutput()


def main():
    """
    Reads the file names from the command options and runs the tests.
    """
    parser = argparse.ArgumentParser(
        description="Compares the contents of two VTU files."
    )
    parser.add_argument("result", type=pathlib.Path, help="the result *.vtu file")
    parser.add_argument("reference", type=pathlib.Path, help="the reference *.vtu file")
    arguments, other = parser.parse_known_args()

    global RESULT_GRID  # pylint: disable=W0603
    RESULT_GRID = read_vtu(arguments.result)

    global REFERENCE_GRID  # pylint: disable=W0603
    REFERENCE_GRID = read_vtu(arguments.reference)

    return unittest.main(argv=sys.argv[:1] + other)


if __name__ == "__main__":
    main()
