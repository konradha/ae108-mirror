#! /usr/bin/env python3

# © 2020 ETH Zurich, Mechanics and Materials Lab
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
Prints quadrature definitions for intervals in 1D and products of intervals in higher dimensions.
"""

import itertools
from typing import Any, Iterable, List, Tuple

import sympy


def compute_gaussian_nodes(order: int) -> Iterable[Any]:
    """
    Computes the nodes of the Gaussian quadrature of the given order.

    >>> compute_gaussian_nodes(1)
    [0]
    >>> compute_gaussian_nodes(2)
    [-sqrt(3)/3, sqrt(3)/3]
    >>> compute_gaussian_nodes(3)
    [-sqrt(15)/5, 0, sqrt(15)/5]
    """

    return sympy.legendre_poly(order, sympy.Symbol("x"), polys=True).real_roots()


def compute_gaussian_weights(order: int) -> List[Any]:
    """
    Computes the weights for the Gaussian quadrature of the given order.

    >>> compute_gaussian_weights(1)
    [2]
    >>> compute_gaussian_weights(2)
    [1, 1]
    >>> compute_gaussian_weights(3)
    [5/9, 8/9, 5/9]
    """

    return [
        2
        / (1 - node**2)
        / sympy.diff(
            sympy.legendre_poly(order, sympy.Symbol("x"), polys=True), sympy.Symbol("x")
        ).eval(node)
        ** 2
        for node in compute_gaussian_nodes(order)
    ]


def compute_gaussian_nodes_in_dimension(
    order: int, dimension: int
) -> List[Tuple[Any, ...]]:
    """
    Computes the nodes of the Gaussian quadrature of the given dimension and order.

    >>> compute_gaussian_nodes_in_dimension(1, 1)
    [(0,)]
    >>> compute_gaussian_nodes_in_dimension(1, 2)
    [(0, 0)]
    >>> compute_gaussian_nodes_in_dimension(2, 1)
    [(-sqrt(3)/3,), (sqrt(3)/3,)]
    >>> len(compute_gaussian_nodes_in_dimension(2, 2))
    4
    >>> compute_gaussian_nodes_in_dimension(2, 2)[0]
    (-sqrt(3)/3, -sqrt(3)/3)
    >>> compute_gaussian_nodes_in_dimension(2, 2)[1]
    (-sqrt(3)/3, sqrt(3)/3)
    >>> compute_gaussian_nodes_in_dimension(2, 2)[2]
    (sqrt(3)/3, -sqrt(3)/3)
    >>> compute_gaussian_nodes_in_dimension(2, 2)[3]
    (sqrt(3)/3, sqrt(3)/3)
    """

    return list(itertools.product(compute_gaussian_nodes(order), repeat=dimension))


def compute_gaussian_weights_in_dimension(order: int, dimension: int) -> List[Any]:
    """
    Computes the nodes of the Gaussian quadrature of the given dimension and order.

    >>> compute_gaussian_weights_in_dimension(1, 1)
    [2.00000000000000]
    >>> compute_gaussian_weights_in_dimension(1, 2)
    [4.00000000000000]
    >>> compute_gaussian_weights_in_dimension(2, 1)
    [1.00000000000000, 1.00000000000000]
    >>> compute_gaussian_weights_in_dimension(2, 2)
    [1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000]
    """

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
    """
    Converts the input to a C++ initializer list.

    >>> to_initializer_list([1., 2.])
    '{{+1.0000000000000000, +2.0000000000000000}}'
    >>> to_initializer_list((1., 2.))
    '{{+1.0000000000000000, +2.0000000000000000}}'
    >>> to_initializer_list(1.)
    '+1.0000000000000000'
    """

    if isinstance(values, (list, tuple)):
        result = ", ".join([to_initializer_list(element) for element in values])
        return "{{" + result + (", " if result.endswith("}") else "") + "}}"
    return f"{float(values):+.16f}"


def quadrature_definition(order: int, dimension: int) -> str:
    """
    Returns the definition of a quadrature rule.

    >>> quadrature_definition(1, 1)
    'AE108_ELEMENTS_QUADRATURE_DEFINE(QuadratureType::Cube, 1, 1, 1, \
{{{{{+0.0000000000000000}}, }}, {{+2.0000000000000000}}});'
    """

    nodes = compute_gaussian_nodes_in_dimension(order, dimension)
    weights = compute_gaussian_weights_in_dimension(order, dimension)
    return (
        "AE108_ELEMENTS_QUADRATURE_DEFINE("
        f"QuadratureType::Cube, {dimension}, {2 * order - 1}, {len(nodes)}, "
        f"{{{to_initializer_list(nodes)}, {to_initializer_list(weights)}}}"
        ");"
    )


def main() -> None:
    """Prints quadrature definitions."""
    for order in range(1, 6):
        for dimension in range(1, 4):
            print(quadrature_definition(order, dimension))
            print()


if __name__ == "__main__":
    main()
