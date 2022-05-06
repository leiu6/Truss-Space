# Truss Space
# 2D FEA Solver for Trusses


from builder import *
import numpy as np


def main():
    # Define material properties
    E = 29.5 * 10 ** 6
    A = 1

    # Initialize structure
    s1 = Structure(E, A)

    # Add elements
    s1.add_element((0, 0), (40, 0))
    s1.add_element((40, 30), (40, 0))
    s1.add_element((0, 0), (40, 30))
    s1.add_element((0, 30), (40, 30))

    # Apply node boundary conditions
    s1.add_boundary_condition(1, (0, 0))
    s1.add_boundary_condition(2, ('free', 0))
    s1.add_boundary_condition(3, ('free', 'free'))
    s1.add_boundary_condition(4, (0, 0))

    # Apply node forces
    s1.apply_force(1, 0, 0)
    s1.apply_force(2, 20000, 0)
    s1.apply_force(3, 0, -25000)
    s1.apply_force(4, 0, 0)

    # Generate global stiffness matrix
    s1.build_global_stiffness_matrix()

    # Solve
    s1.solve_displacements()
    s1.solve_stresses()
    s1.solve_reactions()

    # Display
    print("Displacements:", s1.displacements)
    print("Stresses:", s1.stresses)
    print("Reactions:", s1.reactions)


if __name__ == "__main__":
    main()
