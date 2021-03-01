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
Generates an XDMF file from an ae108 file.
This XDMF file can be opened with Paraview.

The script takes one parameter: the path to the ae108 file.
It writes an XDMF file of the same name to the working directory.
"""

import argparse
from enum import IntEnum
import pathlib
import re
import typing
import sys
import xml.etree.ElementTree as ET

import h5py
import numpy


def parse_input_filename() -> pathlib.Path:
    """
    Parses the command line parameters and returns the input file name.
    """
    parser = argparse.ArgumentParser(
        description="Generate XDMF files for *.ae108 files."
    )
    parser.add_argument("input", help="*.ae108 file to read", type=pathlib.Path)
    return parser.parse_args().input


def extract_field_name(name: str) -> str:
    """
    Removes the trailing "_Field_0" (and similar) from the input name.

    >>> extract_field_name("abc")
    'abc'
    >>> extract_field_name("abc_Field_0")
    'abc'
    >>> extract_field_name("abc_Field_10")
    'abc'
    >>> extract_field_name("abc_field_10")
    'abc_field_10'
    """
    return re.sub(r"_Field_[\d]+$", "", name)


class UnsupportedElementType(Exception):
    """
    This exception is raised if an unsupported element type
    (number of vertices, topological dimension) is provided.
    """


def number_of_corners_to_type(
    number_of_vertices: int, topological_dimension: int
) -> str:
    """
    Returns a guess of the element type for the provided number of vertices.

    >>> number_of_corners_to_type(3, 2)
    'Triangle'
    >>> number_of_corners_to_type(4, 2)
    'Quadrilateral'
    >>> number_of_corners_to_type(6, 2)
    'Tri_6'
    >>> number_of_corners_to_type(2, 2)
    Traceback (most recent call last):
    generate_xdmf.UnsupportedElementType
    >>> number_of_corners_to_type(5, 2)
    Traceback (most recent call last):
    generate_xdmf.UnsupportedElementType
    >>> number_of_corners_to_type(4, 3)
    'Tetrahedron'
    >>> number_of_corners_to_type(8, 3)
    'Hexahedron'
    >>> number_of_corners_to_type(3, 3)
    Traceback (most recent call last):
    generate_xdmf.UnsupportedElementType
    >>> number_of_corners_to_type(9, 3)
    Traceback (most recent call last):
    generate_xdmf.UnsupportedElementType
    >>> number_of_corners_to_type(4, 1)
    Traceback (most recent call last):
    generate_xdmf.UnsupportedElementType
    """
    try:
        type_map = {
            2: {
                3: "Triangle",
                4: "Quadrilateral",
                6: "Tri_6",
            },
            3: {
                4: "Tetrahedron",
                8: "Hexahedron",
            },
        }
        return type_map[topological_dimension][number_of_vertices]
    except KeyError as error:
        raise UnsupportedElementType() from error


class UnsupportedHdfDtype(Exception):
    """
    This exception is raised if an unsupported HDF datatype is provided.
    """


def hdf_dtype_to_type(dtype: numpy.dtype) -> typing.Tuple[str, str, int]:
    """
    Converts a Numpy dtype to a tuple describing the datatype for XDMF.
    The tuple has the following components: endianness, number type, bytes.

    >>> hdf_dtype_to_type(numpy.dtype('<f8'))
    ('Little', 'Float', 8)
    >>> hdf_dtype_to_type(numpy.dtype('<f4'))
    ('Little', 'Float', 4)
    >>> hdf_dtype_to_type(numpy.dtype('<i8'))
    ('Little', 'Int', 8)
    >>> hdf_dtype_to_type(numpy.dtype('<i4'))
    ('Little', 'Int', 4)
    >>> hdf_dtype_to_type(numpy.dtype(numpy.int64))
    ('Little', 'Int', 8)
    >>> hdf_dtype_to_type(numpy.dtype(numpy.uint64))
    Traceback (most recent call last):
    generate_xdmf.UnsupportedHdfDtype
    """
    try:
        datatype_map = {
            numpy.dtype("<f8"): ("Little", "Float", 8),
            numpy.dtype("<f4"): ("Little", "Float", 4),
            numpy.dtype("<i8"): ("Little", "Int", 8),
            numpy.dtype("<i4"): ("Little", "Int", 4),
        }
        return datatype_map[dtype]
    except KeyError as error:
        raise UnsupportedHdfDtype() from error


class UnsupportedCoordinateDimension(Exception):
    """
    This exception is raised if an unsupported coordinate dimension is provided.
    """


def coordinate_dimension_to_type(dimension: int) -> str:
    """
    Returns the grid type for the given dimension.

    >>> coordinate_dimension_to_type(2)
    'XY'
    >>> coordinate_dimension_to_type(3)
    'XYZ'
    >>> coordinate_dimension_to_type(1)
    Traceback (most recent call last):
    generate_xdmf.UnsupportedCoordinateDimension
    >>> coordinate_dimension_to_type(4)
    Traceback (most recent call last):
    generate_xdmf.UnsupportedCoordinateDimension
    """
    try:
        type_map = {
            2: "XY",
            3: "XYZ",
        }
        return type_map[dimension]
    except KeyError as error:
        raise UnsupportedCoordinateDimension() from error


def shape_to_xdmf_dimensions(shape: typing.Tuple[int, int]) -> str:
    """
    Returns the dimension provided by shape

    >>> shape_to_xdmf_dimensions((2, 3))
    '2 3'
    """
    return " ".join(str(x) for x in shape)


def add_hdf_dataitem(
    parent: ET.Element, hdf_file: h5py.File, hdf_path: str
) -> ET.Element:
    """
    Inserts a DataItem child with name that refers to the path `hdf_path`
    in the file at `file_name`.
    """
    item = hdf_file[hdf_path]
    datatype = hdf_dtype_to_type(item.dtype)
    dataitem = ET.SubElement(
        parent,
        "DataItem",
        {
            "Dimensions": shape_to_xdmf_dimensions(hdf_file[hdf_path].shape),
            "Format": "HDF",
            "Endian": str(datatype[0]),
            "NumberType": str(datatype[1]),
            "Precision": str(datatype[2]),
            "Rank": str(item.ndim),
        },
    )
    dataitem.text = "{}:{}".format(hdf_file.filename, hdf_path)
    return dataitem


def read_topological_dimension(hdf_file: h5py.File) -> int:
    """
    Returns the topological dimension specified in the mesh.
    """
    return hdf_file["/viz/topology/cells"].attrs.get("cell_dim")


def read_coordinate_dimension(hdf_file: h5py.File) -> int:
    """
    Returns the topological dimension specified in the mesh.
    """
    return hdf_file["/geometry/vertices"].shape[1]


class MeshInformationMissing(Exception):
    """
    This exception is raised if mesh information is missing in the HDF5 file.
    """


def add_topology(parent: ET.Element, hdf_file: h5py.File) -> ET.Element:
    """
    Adds the topology element to the parent.
    """
    if not "/viz/topology/cells" in hdf_file:
        raise MeshInformationMissing()

    number_of_elements = hdf_file["/viz/topology/cells"].shape[0]
    number_of_corners = hdf_file["/viz/topology/cells"].attrs.get("cell_corners")

    topology = ET.SubElement(
        parent,
        "Topology",
        {
            "NumberOfElements": str(number_of_elements),
            "Type": number_of_corners_to_type(
                number_of_corners, read_topological_dimension(hdf_file)
            ),
        },
    )
    add_hdf_dataitem(topology, hdf_file, "/viz/topology/cells")

    return topology


def add_geometry(parent: ET.Element, hdf_file: h5py.File) -> ET.Element:
    """
    Adds the geometry element to the parent.
    """
    if not "/geometry/vertices" in hdf_file:
        raise MeshInformationMissing()

    coordinate_dimension = read_coordinate_dimension(hdf_file)
    geometry = ET.SubElement(
        parent,
        "Geometry",
        {"GeometryType": coordinate_dimension_to_type(coordinate_dimension)},
    )
    add_hdf_dataitem(geometry, hdf_file, "/geometry/vertices")
    return geometry


def add_column_selector_dataitem(
    parent: ET.Element, hdf_file: h5py.File, hdf_path: str, column: int
) -> ET.Element:
    """
    Adds a dataitem describing selected columns of a two dimensional dataset.
    """
    ndim = hdf_file[hdf_path].ndim
    assert ndim == 2
    shape = hdf_file[hdf_path].shape
    assert column >= 0
    assert column < shape[1]

    dataitem = ET.SubElement(
        parent,
        "DataItem",
        {"Format": "XML", "Dimensions": "3 {}".format(ndim)},
    )
    dataitem.text = "0 {} 1 1 {} 1".format(column, shape[0])
    return dataitem


def add_field(
    parent: ET.Element, hdf_file: h5py.File, field_name: str, center: str
) -> None:
    """
    Adds a field with name `field_name`, centered at `center`, to `parent`.
    If the field only has one dimension, then a scalar field will be added.
    If the field has two dimensions, then there are two cases:
        - 3 columns: A vector field will be added.
        - otherwise: The field is split into scalar fields, and those are added.
    """
    hdf_path = "/fields/{}".format(field_name)

    ndim = hdf_file[hdf_path].ndim
    shape = hdf_file[hdf_path].shape
    assert ndim >= 1
    assert ndim <= 2

    if ndim == 1:
        attribute = ET.SubElement(
            parent,
            "Attribute",
            {
                "Name": field_name,
                "Center": center,
                "AttributeType": "Scalar",
            },
        )
        add_hdf_dataitem(attribute, hdf_file, hdf_path)
        return

    if ndim == 2 and shape[1] == 3:
        attribute = ET.SubElement(
            parent,
            "Attribute",
            {
                "Name": field_name,
                "Center": center,
                "AttributeType": "Vector",
            },
        )
        add_hdf_dataitem(attribute, hdf_file, hdf_path)
        return

    for column in range(0, shape[1]):
        attribute = ET.SubElement(
            parent,
            "Attribute",
            {
                "Name": "{}[{}]".format(field_name, column),
                "Center": center,
                "AttributeType": "Scalar",
            },
        )
        hyperslab = ET.SubElement(
            attribute,
            "DataItem",
            {"ItemType": "HyperSlab", "Dimensions": "{} 1".format(shape[0])},
        )
        add_column_selector_dataitem(hyperslab, hdf_file, hdf_path, column)
        add_hdf_dataitem(hyperslab, hdf_file, hdf_path)


def add_vertex_fields(parent: ET.Element, hdf_file: h5py.File) -> None:
    """
    Adds all the vertex fields to the parent.
    """
    if not "vertex_fields" in hdf_file:
        return

    vertex_field_names = [
        extract_field_name(field_name)
        for field_name in hdf_file["vertex_fields"].keys()
    ]

    for field_name in vertex_field_names:
        add_field(parent, hdf_file, field_name, "Node")


def add_cell_fields(parent: ET.Element, hdf_file: h5py.File) -> None:
    """
    Adds all the cell fields to the parent.
    """
    if not "cell_fields" in hdf_file:
        return

    cell_field_names = [
        extract_field_name(field_name) for field_name in hdf_file["cell_fields"].keys()
    ]

    for field_name in cell_field_names:
        add_field(parent, hdf_file, field_name, "Cell")


def hdf_to_xdmf_string(hdf_file: h5py.File) -> str:
    """
    Returns the XDMF string corresponding to the provided HDF5 file.
    """
    xdmf = ET.Element("Xdmf")
    domain = ET.SubElement(xdmf, "Domain")
    grid = ET.SubElement(domain, "Grid")

    add_topology(grid, hdf_file)
    add_geometry(grid, hdf_file)
    add_vertex_fields(grid, hdf_file)
    add_cell_fields(grid, hdf_file)

    return ET.tostring(xdmf, encoding="unicode")


def main() -> None:
    """
    Parses the input file name from the command line parameters and writes the corresponding
    XDMF to file.
    """
    filename = parse_input_filename()
    with open(filename.stem + ".xdmf", "w") as xdmf_file:
        xdmf_file.write(hdf_to_xdmf_string(h5py.File(filename, "r")))


class ErrorCode(IntEnum):
    """
    The error codes used by this script on exit.
    """

    OK = 0
    MISSING_MESH = 2
    UNSUPPORTED_COORDINATE_DIM = 3
    UNSUPPORTED_ELEMENT_TYPE = 4
    UNSUPPORTED_DATA_TYPE = 5


if __name__ == "__main__":
    try:
        main()
    except MeshInformationMissing:
        print(
            "Error: "
            "Mesh information is missing in the input file. "
            "Please export the mesh (topology, geometry) to this file.",
            file=sys.stderr,
        )
        sys.exit(ErrorCode.MISSING_MESH)
    except UnsupportedCoordinateDimension:
        print(
            "Error: Coordinate dimensions other than 2 and 3 are not supported.",
            file=sys.stderr,
        )
        sys.exit(ErrorCode.UNSUPPORTED_COORDINATE_DIM)
    except UnsupportedElementType:
        print(
            "Error: "
            "The combination of topological dimension and "
            "number of vertices per element is not supported.",
            file=sys.stderr,
        )
        sys.exit(ErrorCode.UNSUPPORTED_ELEMENT_TYPE)
    except UnsupportedHdfDtype:
        print(
            "Error: The input file contains data of an unsupported type.",
            file=sys.stderr,
        )
        sys.exit(ErrorCode.UNSUPPORTED_DATA_TYPE)
