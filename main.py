# Truss Space
# 2D FEA Solver for Trusses


import solver as s
import numpy as np


def main():
    n1 = (0, 0)
    n2 = (40, 0)
    n3 = (40, 30)
    n4 = (0, 30)

    E = 29.5 * 10 ** 6
    A = 1

    k1 = s.line_elem_stiffness(n1, n2, E, A)
    k2 = s.line_elem_stiffness(n3, n2, E, A)
    k3 = s.line_elem_stiffness(n1, n3, E, A)
    k4 = s.line_elem_stiffness(n4, n3, E, A)

    gsm = s.line_elem_global_stiffness_matrix([(1, 2), (3, 2), (1, 3), (4, 3)], [k1, k2, k3, k4])

    F = np.array([
        [0],
        [0],
        [20000],
        [0],
        [0],
        [-25000],
        [0],
        [0]
    ])

    boundary_conditions = [
        (0, 0),
        ('free', 0),
        ('free', 'free'),
        (0, 0)
    ]

    gsm, F = s.line_elem_boundary_conditions(gsm, F, boundary_conditions)

    print(gsm, F)


if __name__ == "__main__":
    main()
