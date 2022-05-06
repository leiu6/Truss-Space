# Truss Space
# 2D FEA Solver for Trusses


import solver as s
import numpy as np


def main():
    # define nodes
    n1 = (0, 0)
    n2 = (40, 0)
    n3 = (40, 30)
    n4 = (0, 30)

    nodes = [
        n1,
        n2,
        n3,
        n4
    ]

    E = 29.5 * 10 ** 6
    A = 1

    # k1 = s.elem_stiffness(n1, n2, E, A)
    # k2 = s.elem_stiffness(n3, n2, E, A)
    # k3 = s.elem_stiffness(n1, n3, E, A)
    # k4 = s.elem_stiffness(n4, n3, E, A)

    starts_ends = [(1, 2), (3, 2), (1, 3), (4, 3)]

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

    kl = s.elem_stiffness_list(starts_ends, nodes, E, A)

    K_unmodified = s.global_stiffness_matrix(starts_ends, kl)
    K = np.copy(K_unmodified)

    K, F = s.apply_boundary_conditions(K, F, boundary_conditions)

    d = s.solve_displacements(K, F)

    print("Displacements:", d)

    sigma = s.solve_stresses(d, nodes, starts_ends, boundary_conditions, E)

    print("Stresses:", sigma)

    reactions = s.solve_reactions(K_unmodified, d, boundary_conditions)

    print("Reactions:", reactions)


if __name__ == "__main__":
    main()
