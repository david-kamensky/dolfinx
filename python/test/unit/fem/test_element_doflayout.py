# Copyright (C) 2019 Chris Richardson
#
# This file is part of DOLFIN (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Unit tests for assembly over domains"""

import numpy

import dolfin
from dolfin import cpp
import ufl


def test_element_dof_layout():
    mesh = dolfin.generation.UnitSquareMesh(dolfin.MPI.comm_world, 4, 4)

    # Simple layout, one dofs per vertex
    P1_layout = cpp.fem.ElementDofLayout(1, [[[0], [1], [2]]],
                                         [[[0], [1], [2]],
                                          [[1, 2], [0, 2], [0, 1]],
                                          [[0, 1, 2]]],
                                         [], [])

    cpp_dofmap = cpp.fem.DofMap(P1_layout, mesh)
    dofmap = dolfin.DofMap(cpp_dofmap)
    for i in range(mesh.num_entities(2)):
        print (dofmap.cell_dofs(i))

    # Vector layout, contiguous
    P1_x = cpp.fem.ElementDofLayout(1, [[[0], [1], [2]]],
                                    [[[0], [1], [2]],
                                     [[1, 2], [0, 2], [0, 1]],
                                     [[0, 1, 2]]],
                                    [0, 1, 2], [])

    P1_y = cpp.fem.ElementDofLayout(1, [[[0], [1], [2]]],
                                    [[[0], [1], [2]],
                                     [[1, 2], [0, 2], [0, 1]],
                                     [[0, 1, 2]]],
                                    [3, 4, 5], [])

    P1_z = cpp.fem.ElementDofLayout(1, [[[0], [1], [2]]],
                                    [[[0], [1], [2]],
                                     [[1, 2], [0, 2], [0, 1]],
                                     [[0, 1, 2]]],
                                    [6, 7, 8], [])


    Vector_layout = cpp.fem.ElementDofLayout(1, [[[0, 3, 6], [1, 4, 7], [2, 5, 8]]],
                                             [[[0, 3, 6], [1, 4, 7], [2, 5, 8]],
                                             [[1, 2, 4, 5, 7, 8],
                                              [0, 2, 3, 5, 6, 8],
                                              [0, 1, 3, 4, 6, 7]],
                                              [[0, 1, 2, 3, 4, 5, 6, 7, 8]]],
                                              [], [P1_x, P1_y, P1_z])

    cpp_dofmap = cpp.fem.DofMap(Vector_layout, mesh)
    dofmap = dolfin.DofMap(cpp_dofmap)

    for i in range(mesh.num_entities(2)):
        print (dofmap.cell_dofs(i))

    # Vector layout, interleaved
    P1_x = cpp.fem.ElementDofLayout(1, [[[0], [1], [2]]],
                                    [[[0], [1], [2]],
                                     [[1, 2], [0, 2], [0, 1]],
                                     [[0, 1, 2]]],
                                    [0, 3, 6], [])

    P1_y = cpp.fem.ElementDofLayout(1, [[[0], [1], [2]]],
                                    [[[0], [1], [2]],
                                     [[1, 2], [0, 2], [0, 1]],
                                     [[0, 1, 2]]],
                                    [1, 4, 7], [])

    P1_z = cpp.fem.ElementDofLayout(1, [[[0], [1], [2]]],
                                    [[[0], [1], [2]],
                                     [[1, 2], [0, 2], [0, 1]],
                                     [[0, 1, 2]]],
                                    [2, 5, 8], [])


    Vector_layout = cpp.fem.ElementDofLayout(1, [[[0, 1, 2], [3, 4, 5], [6, 7, 8]]],
                                             [[[0, 1, 2], [3, 4, 5], [6, 7, 8]],
                                             [[3, 4, 5, 6, 7, 8],
                                              [0, 1, 2, 6, 7, 8],
                                              [0, 1, 2, 3, 4, 5]],
                                              [[0, 1, 2, 3, 4, 5, 6, 7, 8]]],
                                              [], [P1_x, P1_y, P1_z])

    cpp_dofmap = cpp.fem.DofMap(Vector_layout, mesh)
    dofmap = dolfin.DofMap(cpp_dofmap)

    for i in range(mesh.num_entities(2)):
        print (dofmap.cell_dofs(i))
