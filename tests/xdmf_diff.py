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
Compares the contents of two XDMF files.
"""

import unittest
import typing
import pathlib
import argparse
import sys
from functools import partial
import vtk

RESULT_GRID = vtk.vtkUnstructuredGrid()
REFERENCE_GRID = vtk.vtkUnstructuredGrid()


def point_coordinates(
    grid: vtk.vtkUnstructuredGrid, point_id: int
) -> typing.List[float]:
    """
    Returns the coordinates of the cell with index `point_id`.
    """
    return list(grid.GetPoint(point_id))


def cell_connectivity(grid: vtk.vtkUnstructuredGrid, cell_id: int) -> typing.List[int]:
    """
    Returns the connectivity for the cell with index `cell_id`.
    """
    ids = vtk.vtkIdList()
    grid.GetCells().GetCellAtId(cell_id, ids)
    return list(map(ids.GetId, range(ids.GetNumberOfIds())))


def all_point_coordinates(
    grid: vtk.vtkUnstructuredGrid,
) -> typing.List[typing.List[float]]:
    """
    Returns a list of point coordinates.
    """
    return list(map(partial(point_coordinates, grid), range(grid.GetNumberOfPoints())))


def all_cell_coordinates(
    grid: vtk.vtkUnstructuredGrid,
) -> typing.List[typing.List[typing.List[float]]]:
    """
    Returns a list of point coordinates for every element.
    """
    connectivities = (
        cell_connectivity(grid, i) for i in range(grid.GetCells().GetNumberOfCells())
    )

    def id_to_coordinates(i: int) -> typing.List[float]:
        return point_coordinates(grid, i)

    def connectivity_to_list_of_coordinates(
        connectivity: typing.List[int],
    ) -> typing.List[typing.List[float]]:
        return list(map(id_to_coordinates, connectivity))

    return list(map(connectivity_to_list_of_coordinates, connectivities))


def all_point_data(
    grid: vtk.vtkUnstructuredGrid,
) -> typing.List[typing.Dict[str, typing.Any]]:
    """
    Returns a list of point data.
    """
    data = grid.GetPointData()
    return [
        {
            "name": data.GetArrayName(i),
            "array": [
                list(data.GetArray(i).GetTuple(j))
                for j in range(data.GetArray(i).GetNumberOfTuples())
            ],
        }
        for i in range(data.GetNumberOfArrays())
    ]


def all_cell_data(
    grid: vtk.vtkUnstructuredGrid,
) -> typing.List[typing.Dict[str, typing.Any]]:
    """
    Returns a list of cell data.
    """
    data = grid.GetCellData()
    return [
        {
            "name": data.GetArrayName(i),
            "array": [
                list(data.GetArray(i).GetTuple(j))
                for j in range(data.GetArray(i).GetNumberOfTuples())
            ],
        }
        for i in range(data.GetNumberOfArrays())
    ]


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
        for result_coordinates, reference_coordinates in zip(
            sorted(all_point_coordinates(RESULT_GRID)),
            sorted(all_point_coordinates(REFERENCE_GRID)),
        ):
            for result_x, reference_x in zip(result_coordinates, reference_coordinates):
                self.assertAlmostEqual(result_x, reference_x)

    def test_same_cell_coordinates(self):
        """
        The cells of the grids have the vertices at the same coordinates.
        """
        for result_coordinates, reference_coordinates in zip(
            sorted(all_cell_coordinates(RESULT_GRID)),
            sorted(all_cell_coordinates(REFERENCE_GRID)),
        ):
            for result_x, reference_x in zip(result_coordinates, reference_coordinates):
                self.assertAlmostEqual(result_x, reference_x)

    def test_same_number_of_point_data(self):
        """
        The points of the grid have the same number of data.
        """
        self.assertEqual(
            len(all_point_data(RESULT_GRID)), len(all_point_data(REFERENCE_GRID))
        )

    def test_same_point_data_arrays(self):
        """
        The points of the grid have the same point data at the same coordinates.
        """

        for result, reference in zip(
            sorted(all_point_data(RESULT_GRID), key=lambda x: x["name"]),
            sorted(all_point_data(REFERENCE_GRID), key=lambda x: x["name"]),
        ):
            self.assertEqual(result["name"], reference["name"])
            result_array = sorted(
                list(zip(all_point_coordinates(RESULT_GRID), result["array"]))
            )
            reference_array = sorted(
                list(zip(all_point_coordinates(REFERENCE_GRID), reference["array"]))
            )

            for result_item, reference_item in zip(result_array, reference_array):
                for result_x, reference_x in zip(result_item[1], reference_item[1]):
                    self.assertAlmostEqual(result_x, reference_x)

    def test_same_number_of_cell_data(self):
        """
        The cells of the grid have the same number of data.
        """
        self.assertEqual(
            len(all_cell_data(RESULT_GRID)), len(all_cell_data(REFERENCE_GRID))
        )

    def test_same_cell_data_arrays(self):
        """
        The cells of the grid have the same cell data.
        """

        for result, reference in zip(
            sorted(all_cell_data(RESULT_GRID), key=lambda x: x["name"]),
            sorted(all_cell_data(REFERENCE_GRID), key=lambda x: x["name"]),
        ):
            self.assertEqual(result["name"], reference["name"])
            result_array = sorted(
                list(zip(all_cell_coordinates(RESULT_GRID), result["array"]))
            )
            reference_array = sorted(
                list(zip(all_cell_coordinates(REFERENCE_GRID), reference["array"]))
            )

            for result_item, reference_item in zip(result_array, reference_array):
                for result_x, reference_x in zip(result_item[1], reference_item[1]):
                    self.assertAlmostEqual(result_x, reference_x)


def read_xdmf(path: pathlib.Path) -> vtk.vtkUnstructuredGrid:
    """
    Reads the XDMF file at path and returns the result.
    """
    reader = vtk.vtkXdmfReader()
    reader.SetFileName(str(path))
    reader.Update()
    return reader.GetOutput()


def main():
    """
    Reads the file names from the command options and runs the tests.
    """
    parser = argparse.ArgumentParser(
        description="Compares the contents of two XDMF files."
    )
    parser.add_argument("result", type=pathlib.Path, help="the result *.xdmf file")
    parser.add_argument(
        "reference", type=pathlib.Path, help="the reference *.xdmf file"
    )
    arguments, other = parser.parse_known_args()

    global RESULT_GRID  # pylint: disable=W0603
    RESULT_GRID = read_xdmf(arguments.result)

    global REFERENCE_GRID  # pylint: disable=W0603
    REFERENCE_GRID = read_xdmf(arguments.reference)

    return unittest.main(argv=sys.argv[:1] + other)


if __name__ == "__main__":
    main()
