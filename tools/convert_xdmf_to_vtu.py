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
