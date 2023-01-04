#!/usr/bin/env python3

# © 2022 ETH Zurich, Mechanics and Materials Lab
#
# This file is part of ae108.
#
# ae108 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any
# later version.
#
# ae108 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ae108. If not, see <https://www.gnu.org/licenses/>.

"""
Converts an XDMF file to an ASCII VTU file.
"""

import argparse
import pathlib
from vtk import vtkXMLUnstructuredGridWriter, vtkXdmfReader  # pylint: disable=E0611


def parse_input_filename() -> pathlib.Path:
    """
    Parses the command line parameters and returns the input file name.
    """
    parser = argparse.ArgumentParser(
        description="Converts an XDMF file to an ASCII VTU file."
    )
    parser.add_argument("input", help="*.xdmf file to read", type=pathlib.Path)
    return parser.parse_args().input


def main():
    """
    Reads the XDMF file from the path provided via command line argument
    and writes an ASCII VTU file to the same location, but with extension ".vtu".
    """
    in_filename = parse_input_filename()
    out_filename = in_filename.stem + ".vtu"

    reader = vtkXdmfReader()
    reader.SetFileName(in_filename)
    reader.Update()
    writer = vtkXMLUnstructuredGridWriter()
    writer.SetFileName(out_filename)
    writer.SetDataModeToAscii()
    writer.SetInputData(reader.GetOutput())
    writer.Write()


if __name__ == "__main__":
    main()
