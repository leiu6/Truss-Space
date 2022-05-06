# Truss Space
# 2D FEA Solver for Trusses


import solver as s
import numpy as np


class Structure:
    # This class denotes a structure, which is a collection of line elements connected together at nodes. It stores all
    # relevant parameters to use the solver functions to compute stresses, strains, etc.

    nodes = list()
    starts_ends = list()
    boundary_conditions = list()
    element_stiffness_matrices = list()
    gsm = np.empty(shape=(0, 0))
    gsm_unmodified = np.empty(shape=(0, 0))
    F = list()
    F_unmodified = np.empty(shape=(0, 0))

    # Results storage
    displacements = list()
    stresses = list()
    reactions = list()

    def __init__(self, E, A):
        # Constructor for Structure class
        #
        # E: Young's modulus for structure
        # A: cross-sectional area of line elements

        self.E = E
        self.A = A

    def add_element(self, start_node, end_node):
        # Method to add a line element to the structure
        #
        # start_node: the first node of the line element
        # end_node: the second node of the line element

        k = s.elem_stiffness(start_node, end_node, self.E, self.A)
        self.element_stiffness_matrices.append(k)

        num_nodes = len(self.nodes)  # get initial number of nodes

        # We must check if our nodes already exist or not. If it does exist, we don't need to append it to nodes. We
        # must however append it to the starts_stops
        if start_node not in self.nodes:
            self.nodes.append(start_node)
            start_node_num = num_nodes + 1
        else:
            start_node_num = self.nodes.index(start_node) + 1

        if end_node not in self.nodes:
            self.nodes.append(end_node)
            end_node_num = num_nodes + 2
        else:
            end_node_num = self.nodes.index(end_node) + 1

        # append to starts_stops
        self.starts_ends.append((start_node_num, end_node_num))

    def add_boundary_condition(self, node_number, bc):
        # Add a boundary condition to a given node
        #
        # node_number: index of node
        # bc: boundary condition tuple i.e: (0, 0) or ('free', 0), etc.

        self.boundary_conditions.insert(node_number - 1, bc)

    def apply_force(self, node_number, x_force, y_force):
        # Apply a force to a node
        #
        # node_number: index of node
        # x_force: x component of force
        # y_force: y component of force

        x_insert = 2 * (node_number - 1)
        y_insert = 2 * node_number

        if type(self.F) is np.ndarray:
            self.F = self.F.tolist()

        self.F.insert(x_insert, x_force)
        self.F.insert(y_insert, y_force)

    def build_global_stiffness_matrix(self):
        # Construct initial global stiffness matrix

        self.gsm_unmodified = s.global_stiffness_matrix(self.starts_ends, self.element_stiffness_matrices)

        # Create copy to remove 0 dof rows and cols from. Unmodified version needs to be untouched since reaction force
        # calculations rely on it.
        self.gsm = np.copy(self.gsm_unmodified)

        # Convert to Numpy array
        self.F = np.array(self.F).astype(float)

        self.F_unmodified = np.copy(self.F)

        # Apply boundary conditions
        self.gsm, self.F = s.apply_boundary_conditions(self.gsm, self.F, self.boundary_conditions)

    def solve_displacements(self):
        self.displacements = s.solve_displacements(self.gsm, self.F)
        self.displacements = s.displacement_add_zeros(self.displacements, self.boundary_conditions)

        return self.displacements

    def solve_stresses(self):
        self.stresses = s.solve_stresses(self.displacements,
                                         self.nodes,
                                         self.starts_ends,
                                         self.boundary_conditions,
                                         self.E)

        return self.stresses

    def solve_reactions(self):
        self.reactions = s.solve_reactions(self.gsm_unmodified,
                                           self.displacements,
                                           self.boundary_conditions,
                                           self.F_unmodified)

        return self.reactions
