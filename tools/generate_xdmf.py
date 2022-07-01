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
import itertools

import h5py
import numpy


def parse_command_line_arguments() -> typing.Tuple[pathlib.Path, bool]:
    """
    Parses the command line parameters and returns the input file name and whether
    planar elements are preferred.
    """
    parser = argparse.ArgumentParser(
        description="Generate XDMF files for *.ae108 files."
    )
    parser.add_argument("input", help="*.ae108 file to read", type=pathlib.Path)
    parser.add_argument(
        "--planar",
        help="generate planar elements (e.g. quadrilaterals) in ambiguous cases",
        action="store_true",
    )
    result = parser.parse_args()
    return result.input, result.planar


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
    (number of vertices) is provided.
    """


def number_of_corners_to_type(
    number_of_vertices: int, prefer_planar: bool = False
) -> typing.List[int]:
    """
    Returns a guess of the element type for the provided number of vertices.

    Some element types (e.g. polylines) are described by two numbers:
    the type and the number of vertices.

    >>> number_of_corners_to_type(1) # vertex
    [1, 1]
    >>> number_of_corners_to_type(2) # line
    [2, 2]
    >>> number_of_corners_to_type(3) # triangle
    [4]
    >>> number_of_corners_to_type(4) # tetrahedron
    [6]
    >>> number_of_corners_to_type(4, False) # tetrahedron
    [6]
    >>> number_of_corners_to_type(4, True) # quadrilateral
    [5]
    >>> number_of_corners_to_type(6) # quadratic triangle
    [36]
    >>> number_of_corners_to_type(8) # hexahedron
    [9]
    >>> number_of_corners_to_type(8, False) # hexahedron
    [9]
    >>> number_of_corners_to_type(8, True) # quadratic quadrilateral (8)
    [37]
    >>> number_of_corners_to_type(9) # quadratic quadrilateral (9)
    [35]
    >>> number_of_corners_to_type(10) # quadratic tetrahedron (38)
    [38]
    >>> number_of_corners_to_type(20) # quadratic hexahedron (20)
    [48]
    >>> number_of_corners_to_type(24) # quadratic hexahedron (24)
    [49]
    >>> number_of_corners_to_type(27) # quadratic hexahedron (27)
    [50]
    >>> number_of_corners_to_type(5)
    Traceback (most recent call last):
    generate_xdmf.UnsupportedElementType
    """
    try:
        type_map = {
            1: [1, 1],
            2: [2, 2],
            3: [4],
            4: [5 if prefer_planar else 6],
            6: [36],
            8: [37 if prefer_planar else 9],
            9: [35],
            10: [38],
            20: [48],
            24: [49],
            27: [50],
        }
        return type_map[number_of_vertices]
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


def shape_to_xdmf_dimensions(shape: typing.Tuple[int, ...]) -> str:
    """
    Returns the dimension provided by shape

    >>> shape_to_xdmf_dimensions((2, 3))
    '2 3'
    """
    return " ".join(str(x) for x in shape)


def add_sliced_hdf_dataitem(
    parent: ET.Element,
    hdf_file: h5py.File,
    hdf_path: str,
    sliced_dimensions: typing.Dict[int, int],
) -> ET.Element:
    """
    Inserts a HyperSlab child with name that refers to the path `hdf_path`
    in the file at `file_name`.

    The `sliced_dimension` parameter permits to fix dimensions to a selected index.
    Dimensions that are not present in the HDF5 datasets are ignored.
    """
    item = hdf_file[hdf_path]

    if all((key < 0 or key >= item.ndim) for key in sliced_dimensions):
        add_hdf_dataitem(parent, hdf_file, hdf_path)

    hyperslab_shape = tuple(
        item.shape[i] if i not in sliced_dimensions else 1 for i in range(item.ndim)
    )

    hyperslab = ET.SubElement(
        parent,
        "DataItem",
        {
            "ItemType": "HyperSlab",
            "Dimensions": shape_to_xdmf_dimensions(hyperslab_shape),
        },
    )
    dimensions = ET.SubElement(
        hyperslab,
        "DataItem",
        {"Dimensions": shape_to_xdmf_dimensions((3, item.ndim)), "Format": "XML"},
    )
    dimensions.text = f"""{
          shape_to_xdmf_dimensions(
            tuple(sliced_dimensions.get(i, 0) for i in range(item.ndim))
          )
          } {
          shape_to_xdmf_dimensions(tuple(1 for _ in range(item.ndim)))
          } {
          shape_to_xdmf_dimensions(
            tuple(1 if i in sliced_dimensions else item.shape[i] for i in range(item.ndim))
          )
          }"""
    add_hdf_dataitem(hyperslab, hdf_file, hdf_path)
    return hyperslab


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
    dataitem.text = f"{hdf_file.filename}:{hdf_path}"
    return dataitem


def read_topological_dimension(hdf_file: h5py.File) -> int:
    """
    Returns the topological dimension specified in the mesh.
    """
    return hdf_file["/topology/cells"].attrs.get("cell_dim")


def read_coordinate_dimension(hdf_file: h5py.File) -> int:
    """
    Returns the coordinate dimension specified in the mesh.
    """
    return hdf_file["/geometry/vertices"].shape[1]


class MeshInformationMissing(Exception):
    """
    This exception is raised if mesh information is missing in the HDF5 file.
    """


def read_element_vertices(hdf_file: h5py.File) -> typing.Iterator[numpy.ndarray]:
    """
    Returns the vertex indices of the elements.
    """
    if not "/topology" in hdf_file:
        raise MeshInformationMissing()

    cones_data = hdf_file["/topology/cones"][()]
    number_of_elements = numpy.sum(cones_data > 0)
    cells_data = hdf_file["/topology/cells"][()] - number_of_elements

    offset = 0
    for cone_size in filter(lambda x: x > 0, numpy.nditer(cones_data)):
        yield cells_data[offset : offset + cone_size]
        offset += cone_size


def add_topology(
    parent: ET.Element,
    hdf_file: h5py.File,
    prefer_planar: bool,
) -> ET.Element:
    """
    Writes topology data compatible with XDMF to "/topology/mixed"
    in `hdf_file`, and adds the topology element to the parent.
    """
    if not "/topology" in hdf_file:
        raise MeshInformationMissing()

    cones_data = hdf_file["/topology/cones"][()]

    number_of_elements = numpy.sum(cones_data > 0)
    coordinate_dimension = read_coordinate_dimension(hdf_file)

    topology = ET.SubElement(
        parent,
        "Topology",
        {"TopologyType": "Mixed", "NumberOfElements": str(number_of_elements)},
    )

    topology_data = numpy.fromiter(
        itertools.chain.from_iterable(
            (
                itertools.chain(
                    iter(
                        number_of_corners_to_type(
                            vertices.shape[0],
                            prefer_planar or coordinate_dimension <= 2,
                        )
                    ),
                    numpy.nditer(vertices),
                )
                for vertices in read_element_vertices(hdf_file)
            )
        ),
        dtype="int32",
    )

    hdf_file.require_dataset(
        "/topology/mixed", shape=topology_data.shape, dtype=topology_data.dtype
    )
    hdf_file["/topology/mixed"].write_direct(topology_data)

    add_hdf_dataitem(topology, hdf_file, "/topology/mixed")

    return topology


def is_complex(hdf_file: h5py.File, hdf_path: str) -> bool:
    """
    Returns true if and only if the dataset specifies that it is complex.
    """
    return bool(hdf_file[hdf_path].attrs.get("complex", 0))


def add_geometry(parent: ET.Element, hdf_file: h5py.File) -> ET.Element:
    """
    Adds the geometry element to the parent.
    """
    hdf_path = "/geometry/vertices"
    if not hdf_path in hdf_file:
        raise MeshInformationMissing()

    coordinate_dimension = read_coordinate_dimension(hdf_file)
    geometry = ET.SubElement(
        parent,
        "Geometry",
        {"GeometryType": coordinate_dimension_to_type(coordinate_dimension)},
    )
    complex_values = is_complex(hdf_file, hdf_path)
    assert hdf_file[hdf_path].ndim == (3 if complex_values else 2)

    add_sliced_hdf_dataitem(geometry, hdf_file, hdf_path, sliced_dimensions={2: 0})
    return geometry


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
    hdf_path = f"/fields/{field_name}"

    item = hdf_file[hdf_path]
    shape = item.shape
    complex_values = is_complex(hdf_file, hdf_path)

    field_rank = item.ndim - 1 if complex_values else item.ndim

    assert field_rank >= 1
    assert field_rank <= 2

    loop_variables = [(".real", 0), (".imag", 1)] if complex_values else [("", 0)]

    for postfix, index in loop_variables:
        if field_rank == 1:
            attribute = ET.SubElement(
                parent,
                "Attribute",
                {
                    "Name": field_name + postfix,
                    "Center": center,
                    "AttributeType": "Scalar",
                },
            )
            add_sliced_hdf_dataitem(attribute, hdf_file, hdf_path, {1: index})
            continue

        if field_rank == 2 and shape[1] == 3:
            attribute = ET.SubElement(
                parent,
                "Attribute",
                {
                    "Name": field_name + postfix,
                    "Center": center,
                    "AttributeType": "Vector",
                },
            )
            add_sliced_hdf_dataitem(attribute, hdf_file, hdf_path, {2: index})
            continue

        for column in range(0, shape[1]):
            attribute = ET.SubElement(
                parent,
                "Attribute",
                {
                    "Name": f"{field_name}[{column}]{postfix}",
                    "Center": center,
                    "AttributeType": "Scalar",
                },
            )
            add_sliced_hdf_dataitem(
                attribute, hdf_file, hdf_path, {1: column, 2: index}
            )


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


def hdf_to_xdmf_string(hdf_file: h5py.File, prefer_planar: bool) -> str:
    """
    Returns the XDMF string corresponding to the provided HDF5 file.
    """
    xdmf = ET.Element("Xdmf")
    domain = ET.SubElement(xdmf, "Domain")
    grid = ET.SubElement(domain, "Grid")

    add_topology(grid, hdf_file, prefer_planar)
    add_geometry(grid, hdf_file)
    add_vertex_fields(grid, hdf_file)
    add_cell_fields(grid, hdf_file)

    return ET.tostring(xdmf, encoding="unicode")


def main() -> None:
    """
    Parses the input file name from the command line parameters and writes the corresponding
    XDMF to file.
    """
    filename, prefer_planar = parse_command_line_arguments()
    with open(filename.stem + ".xdmf", "w", encoding="utf-8") as xdmf_file:
        xdmf_file.write(
            hdf_to_xdmf_string(
                h5py.File(filename, "r+"),
                prefer_planar,
            )
        )


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
            "Error: The number of vertices per element is not supported.",
            file=sys.stderr,
        )
        sys.exit(ErrorCode.UNSUPPORTED_ELEMENT_TYPE)
    except UnsupportedHdfDtype:
        print(
            "Error: The input file contains data of an unsupported type.",
            file=sys.stderr,
        )
        sys.exit(ErrorCode.UNSUPPORTED_DATA_TYPE)
