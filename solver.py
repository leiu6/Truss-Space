# Truss Space
# 2D FEA Solver for Trusses
#
# All functions for solving trusses
#
# Helpful guide: http://www.unm.edu/~bgreen/ME360/Finite%20Element%20Truss.pdf


import math
import numpy as np


def line_elem_stiffness(start_point, end_point, E, A):
    # calculates the stiffness matrix for a line element
    #
    # start_point: a tuple that defines the first point of the line element
    # end_point: a tuple that defines the second point of the line element
    # E: Young's modulus
    # A: cross-sectional area of element

    x1 = start_point[0]
    y1 = start_point[1]

    x2 = end_point[0]
    y2 = end_point[1]

    L = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)  # line element length
    c = (x2 - x1) / L  # cosine
    s = (y2 - y1) / L  # sine

    # generate Numpy array for the stiffness matrix
    return E * A / L * np.array([
        [c * c, c * s, -c ** 2, -c * s],
        [c * s, s * s, -c * s, -s ** 2],
        [-c ** 2, -c * s, c * c, c * s],
        [-c * s, -s ** 2, c * s, s * s]
    ])


def line_elem_global_stiffness_matrix(nodes, stiffness_matrices):
    # generates a global stiffness matrix from a list of tuples indicating starting and ending nodes, and a
    # list of global stiffness matrices for each of those elements.
    #
    # nodes: a list of tuples for each element: (start node, end node)
    # stiffness_matrices: a list of numpy arrays for the element stiffness matrices

    dof = 2 * len(nodes)  # Max degrees of freedom. Tells the size of the global stiffness matrix

    # initialize it
    gsm = np.zeros((dof, dof))

    for i, k in enumerate(stiffness_matrices):
        # start and end nodes for the element
        start = nodes[i][0]
        end = nodes[i][1]

        # global degrees of freedom
        g_dof = np.array([line_elem_horizontal_dof(start),
                          line_elem_vertical_dof(start),
                          line_elem_horizontal_dof(end),
                          line_elem_vertical_dof(end)])

        # subtract 1 from each so that they correspond to indices
        g_dof = [x - 1 for x in g_dof]

        # create mini 2x2 matrices
        top_left = k[0:2, 0:2]
        top_right = k[0:2, 2:4]
        bottom_left = k[2:4, 0:2]
        bottom_right = k[2:4, 2:4]

        # place them in global matrix
        gsm[g_dof[0]:g_dof[1] + 1, g_dof[0]:g_dof[1] + 1] += top_left
        gsm[g_dof[0]:g_dof[1] + 1, g_dof[2]:g_dof[3] + 1] += top_right
        gsm[g_dof[2]:g_dof[3] + 1, g_dof[0]:g_dof[1] + 1] += bottom_left
        gsm[g_dof[2]:g_dof[3] + 1, g_dof[2]:g_dof[3] + 1] += bottom_right

    return gsm


def line_elem_vertical_dof(node):
    # computes the vertical degree of freedom for a given node of a line element
    #
    # node: the node number of the line element

    return 2 * node


def line_elem_horizontal_dof(node):
    # computes the horizontal degree of freedom for a given node of a line element
    #
    # node: the node number of the line element

    return 2 * node - 1


def line_elem_boundary_conditions(gsm, F, boundary_conditions):
    # removes the rows and columns of the global stiffness in accordance which the nodes that do not have any
    # displacement in a given direction due to their boundary conditions such as fixed supports, rolling supports, etc.
    #
    # gsm: the original global stiffness matrix that we will remove rows from
    # F: column vector of the applied forces to the structure. This will also have elements removed
    # boundary_conditions: a list of tuples in the following format:
    #   (dx, dy)
    #   -   In order from 1 to # of nodes
    #   -   The displacement constraint at that node in the x or y direction
    #       i.e: ('free', 0): the node is constrained to only move in the x direction (rolling support)
    #            (0, 0): the node cannot move in either direction (fixed support)

    for node, bc in enumerate(boundary_conditions):
        # compute degrees of freedom to be deleted
        dof_h = line_elem_horizontal_dof(node + 1)
        dof_v = line_elem_vertical_dof(node + 1)

        # convert to indices
        i_h = dof_h - 1
        i_v = dof_v - 1

        if bc[0] == 0:
            # horizontal delete
            gsm[i_h, :] = np.nan * np.ones(gsm[i_h, :].size)

            # horizontal delete for F
            F[i_h] = np.nan * np.ones((1, 1))

            # vertical delete
            gsm[:, i_h] = np.nan * np.ones(gsm[:, i_h].size)

        if bc[1] == 0:
            # horizontal delete
            gsm[i_v, :] = np.nan * np.ones(gsm[i_v, :].size)

            # horizontal delete for F
            F[i_v] = np.nan * np.ones((1, 1))

            # vertical delete
            gsm[:, i_v] = np.nan * np.ones(gsm[:, i_v].size)

    # now we delete any nan entries
    gsm = gsm[:, ~np.isnan(gsm).all(axis=0)]
    gsm = gsm[~np.isnan(gsm).all(axis=1), :]
    F = F[F != -2147483648]  # janky workaround. FIX LATER! xD

    return gsm, F
