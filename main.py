# Truss Space
# 2D FEA Solver for Trusses


import solver as s


def main():
    print(s.line_elem_stiffness((1, 1), (2, 2), 20, 1))


if __name__ == "__main__":
    main()
