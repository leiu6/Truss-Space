# Truss Space
# 2D FEA Solver for Trusses
#
# All functions for solving matrices


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
    y2 = end_point[0]

    L = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)  # line element length
    c = (x2 - y1) / L  # cosine
    s = (y2 - y1) / L  # sine

    # generate Numpy array for the stiffness matrix
    return E * A / L * np.array([
        [c * c, c * s, -c ** 2, -c * s],
        [c * s, s * s, -c * s, -s ** 2],
        [-c ** 2, -c * s, c * c, c * s],
        [-c * s, -s ** 2, c * s, s * s]
    ])
