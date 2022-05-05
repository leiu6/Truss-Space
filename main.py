# Truss Space
# 2D FEA Solver for Trusses


import solver as s


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

    print(gsm)


if __name__ == "__main__":
    main()
