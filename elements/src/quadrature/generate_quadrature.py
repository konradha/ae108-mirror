#! /usr/bin/env python3

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
Prints quadrature definitions for intervals in 1D and products of intervals in higher dimensions.
"""

import itertools
from typing import Any, Iterable, List, Tuple

import sympy


def compute_gaussian_nodes(order: int) -> Iterable[Any]:
    """Computes the nodes of the Gaussian quadrature of the given order."""

    return sympy.legendre_poly(order, sympy.Symbol("x"), polys=True).real_roots()


def compute_gaussian_weights(order: int) -> List[Any]:
    """Computes the weights for the Gaussian quadrature of the given order."""

    return [
        2
        / (1 - node ** 2)
        / sympy.diff(
            sympy.legendre_poly(order, sympy.Symbol("x"), polys=True), sympy.Symbol("x")
        ).eval(node)
        ** 2
        for node in compute_gaussian_nodes(order)
    ]


def compute_gaussian_nodes_in_dimension(
    order: int, dimension: int
) -> List[Tuple[Any, ...]]:
    """Computes the nodes of the Gaussian quadrature of the given dimension and order."""

    return list(itertools.product(compute_gaussian_nodes(order), repeat=dimension))


def compute_gaussian_weights_in_dimension(order: int, dimension: int) -> List[Any]:
    """Computes the nodes of the Gaussian quadrature of the given dimension and order."""

    def tuple_product(input_tuple: Tuple[Any, ...]) -> Any:
        product = 1.0
        for value in input_tuple:
            product *= value
        return product

    return [
        tuple_product(product_tuple)
        for product_tuple in itertools.product(
            compute_gaussian_weights(order), repeat=dimension
        )
    ]


def to_initializer_list(values: Any) -> str:
    """Converts the input to a C++ initializer list."""

    if isinstance(values, (list, tuple)):
        result = ", ".join([to_initializer_list(element) for element in values])
        return "{{" + result + (", " if result.endswith("}") else "") + "}}"
    return "{:+.16f}".format(float(values))


def quadrature_definition(order: int, dimension: int) -> str:
    """Returns the definition of a quadrature rule."""

    nodes = compute_gaussian_nodes_in_dimension(order, dimension)
    weights = compute_gaussian_weights_in_dimension(order, dimension)
    return (
        "AE108_ELEMENTS_QUADRATURE_DEFINE(QuadratureType::Cube, "
        + f"{dimension}, {order}, {len(nodes)}, "
        + f"{{{to_initializer_list(nodes)}, {to_initializer_list(weights)}}});"
    )


def main() -> None:
    """Prints quadrature definitions."""
    for order in range(1, 6):
        for dimension in range(1, 4):
            print(quadrature_definition(order, dimension))
            print()


if __name__ == "__main__":
    main()
