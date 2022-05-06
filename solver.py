# Truss Space
# 2D FEA Solver for Trusses
#
# All functions for solving trusses
#
# Helpful guide: http://www.unm.edu/~bgreen/ME360/Finite%20Element%20Truss.pdf


import math
import numpy as np
import scipy.linalg as linalg


def elem_stiffness(start_point, end_point, E, A):
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


def elem_stiffness_list(starts_ends, nodes, E, A):
    # Creates a list of element stiffness matrices for every single element in the structure
    #
    # starts_ends: list of tuples in the format (start_node, end_node)
    # nodes: list of tuples with the initial coordinates for each node
    # E: Young's Modulus
    # A: Cross-sectional area of line elements

    k = list()

    for i, start_end in enumerate(starts_ends):
        start = start_end[0]
        end = start_end[1]

        node_start = nodes[start - 1]
        node_end = nodes[end - 1]

        k.append(elem_stiffness(node_start, node_end, E, A))

    return k


def global_stiffness_matrix(starts_ends, stiffness_matrices):
    # generates a global stiffness matrix from a list of tuples indicating starting and ending nodes, and a
    # list of global stiffness matrices for each of those elements.
    #
    # starts_ends: a list of tuples for each element: (start node, end node)
    # stiffness_matrices: a list of numpy arrays for the element stiffness matrices

    dof = 2 * len(starts_ends)  # Max degrees of freedom. Tells the size of the global stiffness matrix

    # initialize it
    gsm = np.zeros((dof, dof))

    for i, k in enumerate(stiffness_matrices):
        # start and end nodes for the element
        start = starts_ends[i][0]
        end = starts_ends[i][1]

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


def apply_boundary_conditions(gsm, F, boundary_conditions):
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
            F[i_h] = np.nan

            # vertical delete
            gsm[:, i_h] = np.nan * np.ones(gsm[:, i_h].size)

        if bc[1] == 0:
            # horizontal delete
            gsm[i_v, :] = np.nan * np.ones(gsm[i_v, :].size)

            # horizontal delete for F
            F[i_v] = np.nan

            # vertical delete
            gsm[:, i_v] = np.nan * np.ones(gsm[:, i_v].size)

    # now we delete any nan entries
    gsm = gsm[:, ~np.isnan(gsm).all(axis=0)]
    gsm = gsm[~np.isnan(gsm).all(axis=1), :]
    F = F[~np.isnan(F)]

    return gsm, F


def solve_displacements(gsm, F):
    # this function computes the nodal displacements for the entire system
    #
    # gsm: global stiffness matrix after adjustment by removing zero displacement rows and columns
    # F: applied force column vector after removing entries for zero displacement

    # compute displacements
    d = linalg.solve(gsm, F)

    return d


def solve_stresses(displacement, nodes, starts_ends, boundary_conditions, E):
    # computes the normal stress in each line element in the system
    #
    # displacement: vector of the nodal displacements
    # nodes: a list of points for each node
    # starts_ends: a list of tuples noting the start and end nodes of each element. This is again output by the element
    #              stiffness matrix function.
    # E: Young's modulus

    stresses = list()

    trig_vectors, lengths = trig_vectors_and_lengths(nodes, starts_ends)

    for i, trig_vector in enumerate(trig_vectors):
        L = lengths[i]
        d = np.copy(displacement)

        start = starts_ends[i][0]
        end = starts_ends[i][1]

        g_dof = np.array([line_elem_horizontal_dof(start) - 1,
                          line_elem_vertical_dof(start) - 1,
                          line_elem_horizontal_dof(end) - 1,
                          line_elem_vertical_dof(end) - 1])

        # reduce to just the degrees of freedom relevant to given element
        d_trunc = [d[g_dof[0]], d[g_dof[1]], d[g_dof[2]], d[g_dof[3]]]

        # calculate the elemental stress
        sigma = E / L * np.dot(trig_vector, np.transpose(d_trunc))

        stresses.append(sigma)

    return np.array(stresses)


def trig_vectors_and_lengths(nodes, starts_ends):
    # Helper function to assist in calculating stresses. Calculates the trig vectors that are used in the dot product
    # operation with the displacements and also returns the length of the element.
    #
    # start_point: tuple showing the start (x, y) of the line_element
    # end_point: tuple showing the end (x, y) of the line_element

    lengths = list()
    trig_vectors = list()

    for start_end in starts_ends:
        start_node = start_end[0]
        end_node = start_end[1]

        x1 = nodes[start_node - 1][0]
        y1 = nodes[start_node - 1][1]

        x2 = nodes[end_node - 1][0]
        y2 = nodes[end_node - 1][1]

        L = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)  # line element length
        c = (x2 - x1) / L  # cosine
        s = (y2 - y1) / L  # sine

        lengths.append(L)
        trig_vectors.append(np.array([-c, -s, c, s]))

    return trig_vectors, lengths


def displacement_add_zeros(d, boundary_conditions):
    # Add zeros around the computed values to account for the removed rows during the boundary condition application
    # process.
    #
    # d: the displacements output
    # boundary_conditions: the displacement boundary conditions at each node

    new_d = []  # initialize the new list
    j = 0  # number of times we have appended from d

    for i, bc in enumerate(boundary_conditions):
        if bc[0] == 'free':
            new_d.append(d[j])
            j += 1
        else:
            new_d.append(0)

        if bc[1] == 'free':
            new_d.append(d[j])
            j += 1
        else:
            new_d.append(0)

    return np.array(new_d)


def solve_reactions(gsm, d, boundary_conditions, F):
    # Solve for the reaction forces in the structure
    #
    # gsm: global stiffness matrix. Unmodified with no boundary conditions applied.
    # d: displacements vector
    # F: forces vector

    num_nodes = len(d) / 2  # calculate total number of nodes

    reactions = np.matmul(gsm, d)

    # now we must remove entries that are not free elements

    # first we will put all tuple entries in order into a list
    bc_list = list()

    for bc in boundary_conditions:
        bc_list.append(bc[0])
        bc_list.append(bc[1])

    # for i, bc in enumerate(bc_list):
    #     if bc == 'free':
    #         reactions[i] = 0

    return np.subtract(reactions, F)
