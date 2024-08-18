#!/usr/bin/env python3

"""Provides functions for performing the 3D Marching Cubes algorithm.

The Marching Cubes algorithm is a widely used technique for extracting a polygonal mesh 
of an isosurface from a three-dimensional scalar field. This script includes functions 
to perform the Marching Cubes algorithm on a grid of cells and generate meshes representing 
the isosurface within a specified 3D region. It also includes utility functions to create 
obj files for visualization.

Importing necessary modules and defining constants:
"""

from common import adapt, frange
from settings import XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, CELL_SIZE
import math
from utils_3d import V3, Tri, Mesh, make_obj

# Define the vertices of a unit cube in 3D space.
# These vertices are used to define the corners of each cell in the grid.
VERTICES = [
    (0, 0, 0),  # Vertex 0
    (1, 0, 0),  # Vertex 1
    (1, 1, 0),  # Vertex 2
    (0, 1, 0),  # Vertex 3
    (0, 0, 1),  # Vertex 4
    (1, 0, 1),  # Vertex 5
    (1, 1, 1),  # Vertex 6
    (0, 1, 1),  # Vertex 7
]

# Define the edges of the unit cube by specifying pairs of vertices.
# These edges are used to interpolate between vertices to find the boundary of the isosurface.
EDGES = [
    (0, 1),  # Edge 0
    (1, 2),  # Edge 1
    (2, 3),  # Edge 2
    (3, 0),  # Edge 3
    (4, 5),  # Edge 4
    (5, 6),  # Edge 5
    (6, 7),  # Edge 6
    (7, 4),  # Edge 7
    (0, 4),  # Edge 8
    (1, 5),  # Edge 9
    (2, 6),  # Edge 10
    (3, 7),  # Edge 11
]

# Table driven approach to the 256 combinations. Pro-tip, don't write this by hand, copy mine!
# See marching_cubes_gen.py for how I generated these.
# Each index is the bitwise representation of what is solid.
# Each value is a list of triples indicating what edges are used for that triangle
# (Recall each edge of the cell may become a vertex in the output boundary)
cases = [[],
 [[8, 0, 3]],
 [[1, 0, 9]],
 [[8, 1, 3], [8, 9, 1]],
 [[10, 2, 1]],
 [[8, 0, 3], [1, 10, 2]],
 [[9, 2, 0], [9, 10, 2]],
 [[3, 8, 2], [2, 8, 10], [10, 8, 9]],
 [[3, 2, 11]],
 [[0, 2, 8], [2, 11, 8]],
 [[1, 0, 9], [2, 11, 3]],
 [[2, 9, 1], [11, 9, 2], [8, 9, 11]],
 [[3, 10, 11], [3, 1, 10]],
 [[1, 10, 0], [0, 10, 8], [8, 10, 11]],
 [[0, 11, 3], [9, 11, 0], [10, 11, 9]],
 [[8, 9, 11], [11, 9, 10]],
 [[7, 4, 8]],
 [[3, 7, 0], [7, 4, 0]],
 [[7, 4, 8], [9, 1, 0]],
 [[9, 1, 4], [4, 1, 7], [7, 1, 3]],
 [[7, 4, 8], [2, 1, 10]],
 [[4, 3, 7], [4, 0, 3], [2, 1, 10]],
 [[2, 0, 10], [0, 9, 10], [7, 4, 8]],
 [[9, 10, 4], [4, 10, 3], [3, 10, 2], [4, 3, 7]],
 [[4, 8, 7], [3, 2, 11]],
 [[7, 4, 11], [11, 4, 2], [2, 4, 0]],
 [[1, 0, 9], [2, 11, 3], [8, 7, 4]],
 [[2, 11, 1], [1, 11, 9], [9, 11, 7], [9, 7, 4]],
 [[10, 11, 1], [11, 3, 1], [4, 8, 7]],
 [[4, 0, 7], [7, 0, 10], [0, 1, 10], [7, 10, 11]],
 [[7, 4, 8], [0, 11, 3], [9, 11, 0], [10, 11, 9]],
 [[4, 11, 7], [9, 11, 4], [10, 11, 9]],
 [[9, 4, 5]],
 [[9, 4, 5], [0, 3, 8]],
 [[0, 5, 1], [0, 4, 5]],
 [[4, 3, 8], [5, 3, 4], [1, 3, 5]],
 [[5, 9, 4], [10, 2, 1]],
 [[8, 0, 3], [1, 10, 2], [4, 5, 9]],
 [[10, 4, 5], [2, 4, 10], [0, 4, 2]],
 [[3, 10, 2], [8, 10, 3], [5, 10, 8], [4, 5, 8]],
 [[9, 4, 5], [11, 3, 2]],
 [[11, 0, 2], [11, 8, 0], [9, 4, 5]],
 [[5, 1, 4], [1, 0, 4], [11, 3, 2]],
 [[5, 1, 4], [4, 1, 11], [1, 2, 11], [4, 11, 8]],
 [[3, 10, 11], [3, 1, 10], [5, 9, 4]],
 [[9, 4, 5], [1, 10, 0], [0, 10, 8], [8, 10, 11]],
 [[5, 0, 4], [11, 0, 5], [11, 3, 0], [10, 11, 5]],
 [[5, 10, 4], [4, 10, 8], [8, 10, 11]],
 [[9, 7, 5], [9, 8, 7]],
 [[0, 5, 9], [3, 5, 0], [7, 5, 3]],
 [[8, 7, 0], [0, 7, 1], [1, 7, 5]],
 [[7, 5, 3], [3, 5, 1]],
 [[7, 5, 8], [5, 9, 8], [2, 1, 10]],
 [[10, 2, 1], [0, 5, 9], [3, 5, 0], [7, 5, 3]],
 [[8, 2, 0], [5, 2, 8], [10, 2, 5], [7, 5, 8]],
 [[2, 3, 10], [10, 3, 5], [5, 3, 7]],
 [[9, 7, 5], [9, 8, 7], [11, 3, 2]],
 [[0, 2, 9], [9, 2, 7], [7, 2, 11], [9, 7, 5]],
 [[3, 2, 11], [8, 7, 0], [0, 7, 1], [1, 7, 5]],
 [[11, 1, 2], [7, 1, 11], [5, 1, 7]],
 [[3, 1, 11], [11, 1, 10], [8, 7, 9], [9, 7, 5]],
 [[11, 7, 0], [7, 5, 0], [5, 9, 0], [10, 11, 0], [1, 10, 0]],
 [[0, 5, 10], [0, 7, 5], [0, 8, 7], [0, 10, 11], [0, 11, 3]],
 [[10, 11, 5], [11, 7, 5]],
 [[5, 6, 10]],
 [[8, 0, 3], [10, 5, 6]],
 [[0, 9, 1], [5, 6, 10]],
 [[8, 1, 3], [8, 9, 1], [10, 5, 6]],
 [[1, 6, 2], [1, 5, 6]],
 [[6, 2, 5], [2, 1, 5], [8, 0, 3]],
 [[5, 6, 9], [9, 6, 0], [0, 6, 2]],
 [[5, 8, 9], [2, 8, 5], [3, 8, 2], [6, 2, 5]],
 [[3, 2, 11], [10, 5, 6]],
 [[0, 2, 8], [2, 11, 8], [5, 6, 10]],
 [[3, 2, 11], [0, 9, 1], [10, 5, 6]],
 [[5, 6, 10], [2, 9, 1], [11, 9, 2], [8, 9, 11]],
 [[11, 3, 6], [6, 3, 5], [5, 3, 1]],
 [[11, 8, 6], [6, 8, 1], [1, 8, 0], [6, 1, 5]],
 [[5, 0, 9], [6, 0, 5], [3, 0, 6], [11, 3, 6]],
 [[6, 9, 5], [11, 9, 6], [8, 9, 11]],
 [[7, 4, 8], [6, 10, 5]],
 [[3, 7, 0], [7, 4, 0], [10, 5, 6]],
 [[7, 4, 8], [6, 10, 5], [9, 1, 0]],
 [[5, 6, 10], [9, 1, 4], [4, 1, 7], [7, 1, 3]],
 [[1, 6, 2], [1, 5, 6], [7, 4, 8]],
 [[6, 1, 5], [2, 1, 6], [0, 7, 4], [3, 7, 0]],
 [[4, 8, 7], [5, 6, 9], [9, 6, 0], [0, 6, 2]],
 [[2, 3, 9], [3, 7, 9], [7, 4, 9], [6, 2, 9], [5, 6, 9]],
 [[2, 11, 3], [7, 4, 8], [10, 5, 6]],
 [[6, 10, 5], [7, 4, 11], [11, 4, 2], [2, 4, 0]],
 [[1, 0, 9], [8, 7, 4], [3, 2, 11], [5, 6, 10]],
 [[1, 2, 9], [9, 2, 11], [9, 11, 4], [4, 11, 7], [5, 6, 10]],
 [[7, 4, 8], [11, 3, 6], [6, 3, 5], [5, 3, 1]],
 [[11, 0, 1], [11, 4, 0], [11, 7, 4], [11, 1, 5], [11, 5, 6]],
 [[6, 9, 5], [0, 9, 6], [11, 0, 6], [3, 0, 11], [4, 8, 7]],
 [[5, 6, 9], [9, 6, 11], [9, 11, 7], [9, 7, 4]],
 [[4, 10, 9], [4, 6, 10]],
 [[10, 4, 6], [10, 9, 4], [8, 0, 3]],
 [[1, 0, 10], [10, 0, 6], [6, 0, 4]],
 [[8, 1, 3], [6, 1, 8], [6, 10, 1], [4, 6, 8]],
 [[9, 2, 1], [4, 2, 9], [6, 2, 4]],
 [[3, 8, 0], [9, 2, 1], [4, 2, 9], [6, 2, 4]],
 [[0, 4, 2], [2, 4, 6]],
 [[8, 2, 3], [4, 2, 8], [6, 2, 4]],
 [[4, 10, 9], [4, 6, 10], [2, 11, 3]],
 [[11, 8, 2], [2, 8, 0], [6, 10, 4], [4, 10, 9]],
 [[2, 11, 3], [1, 0, 10], [10, 0, 6], [6, 0, 4]],
 [[8, 4, 1], [4, 6, 1], [6, 10, 1], [11, 8, 1], [2, 11, 1]],
 [[3, 1, 11], [11, 1, 4], [1, 9, 4], [11, 4, 6]],
 [[6, 11, 1], [11, 8, 1], [8, 0, 1], [4, 6, 1], [9, 4, 1]],
 [[3, 0, 11], [11, 0, 6], [6, 0, 4]],
 [[4, 11, 8], [4, 6, 11]],
 [[6, 8, 7], [10, 8, 6], [9, 8, 10]],
 [[3, 7, 0], [0, 7, 10], [7, 6, 10], [0, 10, 9]],
 [[1, 6, 10], [0, 6, 1], [7, 6, 0], [8, 7, 0]],
 [[10, 1, 6], [6, 1, 7], [7, 1, 3]],
 [[9, 8, 1], [1, 8, 6], [6, 8, 7], [1, 6, 2]],
 [[9, 7, 6], [9, 3, 7], [9, 0, 3], [9, 6, 2], [9, 2, 1]],
 [[7, 6, 8], [8, 6, 0], [0, 6, 2]],
 [[3, 6, 2], [3, 7, 6]],
 [[3, 2, 11], [6, 8, 7], [10, 8, 6], [9, 8, 10]],
 [[7, 9, 0], [7, 10, 9], [7, 6, 10], [7, 0, 2], [7, 2, 11]],
 [[0, 10, 1], [6, 10, 0], [8, 6, 0], [7, 6, 8], [2, 11, 3]],
 [[1, 6, 10], [7, 6, 1], [11, 7, 1], [2, 11, 1]],
 [[1, 9, 6], [9, 8, 6], [8, 7, 6], [3, 1, 6], [11, 3, 6]],
 [[9, 0, 1], [11, 7, 6]],
 [[0, 11, 3], [6, 11, 0], [7, 6, 0], [8, 7, 0]],
 [[7, 6, 11]],
 [[11, 6, 7]],
 [[3, 8, 0], [11, 6, 7]],
 [[1, 0, 9], [6, 7, 11]],
 [[1, 3, 9], [3, 8, 9], [6, 7, 11]],
 [[10, 2, 1], [6, 7, 11]],
 [[10, 2, 1], [3, 8, 0], [6, 7, 11]],
 [[9, 2, 0], [9, 10, 2], [11, 6, 7]],
 [[11, 6, 7], [3, 8, 2], [2, 8, 10], [10, 8, 9]],
 [[2, 6, 3], [6, 7, 3]],
 [[8, 6, 7], [0, 6, 8], [2, 6, 0]],
 [[7, 2, 6], [7, 3, 2], [1, 0, 9]],
 [[8, 9, 7], [7, 9, 2], [2, 9, 1], [7, 2, 6]],
 [[6, 1, 10], [7, 1, 6], [3, 1, 7]],
 [[8, 0, 7], [7, 0, 6], [6, 0, 1], [6, 1, 10]],
 [[7, 3, 6], [6, 3, 9], [3, 0, 9], [6, 9, 10]],
 [[7, 8, 6], [6, 8, 10], [10, 8, 9]],
 [[8, 11, 4], [11, 6, 4]],
 [[11, 0, 3], [6, 0, 11], [4, 0, 6]],
 [[6, 4, 11], [4, 8, 11], [1, 0, 9]],
 [[1, 3, 9], [9, 3, 6], [3, 11, 6], [9, 6, 4]],
 [[8, 11, 4], [11, 6, 4], [1, 10, 2]],
 [[1, 10, 2], [11, 0, 3], [6, 0, 11], [4, 0, 6]],
 [[2, 9, 10], [0, 9, 2], [4, 11, 6], [8, 11, 4]],
 [[3, 4, 9], [3, 6, 4], [3, 11, 6], [3, 9, 10], [3, 10, 2]],
 [[3, 2, 8], [8, 2, 4], [4, 2, 6]],
 [[2, 4, 0], [6, 4, 2]],
 [[0, 9, 1], [3, 2, 8], [8, 2, 4], [4, 2, 6]],
 [[1, 2, 9], [9, 2, 4], [4, 2, 6]],
 [[10, 3, 1], [4, 3, 10], [4, 8, 3], [6, 4, 10]],
 [[10, 0, 1], [6, 0, 10], [4, 0, 6]],
 [[3, 10, 6], [3, 9, 10], [3, 0, 9], [3, 6, 4], [3, 4, 8]],
 [[9, 10, 4], [10, 6, 4]],
 [[9, 4, 5], [7, 11, 6]],
 [[9, 4, 5], [7, 11, 6], [0, 3, 8]],
 [[0, 5, 1], [0, 4, 5], [6, 7, 11]],
 [[11, 6, 7], [4, 3, 8], [5, 3, 4], [1, 3, 5]],
 [[1, 10, 2], [9, 4, 5], [6, 7, 11]],
 [[8, 0, 3], [4, 5, 9], [10, 2, 1], [11, 6, 7]],
 [[7, 11, 6], [10, 4, 5], [2, 4, 10], [0, 4, 2]],
 [[8, 2, 3], [10, 2, 8], [4, 10, 8], [5, 10, 4], [11, 6, 7]],
 [[2, 6, 3], [6, 7, 3], [9, 4, 5]],
 [[5, 9, 4], [8, 6, 7], [0, 6, 8], [2, 6, 0]],
 [[7, 3, 6], [6, 3, 2], [4, 5, 0], [0, 5, 1]],
 [[8, 1, 2], [8, 5, 1], [8, 4, 5], [8, 2, 6], [8, 6, 7]],
 [[9, 4, 5], [6, 1, 10], [7, 1, 6], [3, 1, 7]],
 [[7, 8, 6], [6, 8, 0], [6, 0, 10], [10, 0, 1], [5, 9, 4]],
 [[3, 0, 10], [0, 4, 10], [4, 5, 10], [7, 3, 10], [6, 7, 10]],
 [[8, 6, 7], [10, 6, 8], [5, 10, 8], [4, 5, 8]],
 [[5, 9, 6], [6, 9, 11], [11, 9, 8]],
 [[11, 6, 3], [3, 6, 0], [0, 6, 5], [0, 5, 9]],
 [[8, 11, 0], [0, 11, 5], [5, 11, 6], [0, 5, 1]],
 [[6, 3, 11], [5, 3, 6], [1, 3, 5]],
 [[10, 2, 1], [5, 9, 6], [6, 9, 11], [11, 9, 8]],
 [[3, 11, 0], [0, 11, 6], [0, 6, 9], [9, 6, 5], [1, 10, 2]],
 [[0, 8, 5], [8, 11, 5], [11, 6, 5], [2, 0, 5], [10, 2, 5]],
 [[11, 6, 3], [3, 6, 5], [3, 5, 10], [3, 10, 2]],
 [[3, 9, 8], [6, 9, 3], [5, 9, 6], [2, 6, 3]],
 [[9, 6, 5], [0, 6, 9], [2, 6, 0]],
 [[6, 5, 8], [5, 1, 8], [1, 0, 8], [2, 6, 8], [3, 2, 8]],
 [[2, 6, 1], [6, 5, 1]],
 [[6, 8, 3], [6, 9, 8], [6, 5, 9], [6, 3, 1], [6, 1, 10]],
 [[1, 10, 0], [0, 10, 6], [0, 6, 5], [0, 5, 9]],
 [[3, 0, 8], [6, 5, 10]],
 [[10, 6, 5]],
 [[5, 11, 10], [5, 7, 11]],
 [[5, 11, 10], [5, 7, 11], [3, 8, 0]],
 [[11, 10, 7], [10, 5, 7], [0, 9, 1]],
 [[5, 7, 10], [10, 7, 11], [9, 1, 8], [8, 1, 3]],
 [[2, 1, 11], [11, 1, 7], [7, 1, 5]],
 [[3, 8, 0], [2, 1, 11], [11, 1, 7], [7, 1, 5]],
 [[2, 0, 11], [11, 0, 5], [5, 0, 9], [11, 5, 7]],
 [[2, 9, 5], [2, 8, 9], [2, 3, 8], [2, 5, 7], [2, 7, 11]],
 [[10, 3, 2], [5, 3, 10], [7, 3, 5]],
 [[10, 0, 2], [7, 0, 10], [8, 0, 7], [5, 7, 10]],
 [[0, 9, 1], [10, 3, 2], [5, 3, 10], [7, 3, 5]],
 [[7, 8, 2], [8, 9, 2], [9, 1, 2], [5, 7, 2], [10, 5, 2]],
 [[3, 1, 7], [7, 1, 5]],
 [[0, 7, 8], [1, 7, 0], [5, 7, 1]],
 [[9, 5, 0], [0, 5, 3], [3, 5, 7]],
 [[5, 7, 9], [7, 8, 9]],
 [[4, 10, 5], [8, 10, 4], [11, 10, 8]],
 [[3, 4, 0], [10, 4, 3], [10, 5, 4], [11, 10, 3]],
 [[1, 0, 9], [4, 10, 5], [8, 10, 4], [11, 10, 8]],
 [[4, 3, 11], [4, 1, 3], [4, 9, 1], [4, 11, 10], [4, 10, 5]],
 [[1, 5, 2], [2, 5, 8], [5, 4, 8], [2, 8, 11]],
 [[5, 4, 11], [4, 0, 11], [0, 3, 11], [1, 5, 11], [2, 1, 11]],
 [[5, 11, 2], [5, 8, 11], [5, 4, 8], [5, 2, 0], [5, 0, 9]],
 [[5, 4, 9], [2, 3, 11]],
 [[3, 4, 8], [2, 4, 3], [5, 4, 2], [10, 5, 2]],
 [[5, 4, 10], [10, 4, 2], [2, 4, 0]],
 [[2, 8, 3], [4, 8, 2], [10, 4, 2], [5, 4, 10], [0, 9, 1]],
 [[4, 10, 5], [2, 10, 4], [1, 2, 4], [9, 1, 4]],
 [[8, 3, 4], [4, 3, 5], [5, 3, 1]],
 [[1, 5, 0], [5, 4, 0]],
 [[5, 0, 9], [3, 0, 5], [8, 3, 5], [4, 8, 5]],
 [[5, 4, 9]],
 [[7, 11, 4], [4, 11, 9], [9, 11, 10]],
 [[8, 0, 3], [7, 11, 4], [4, 11, 9], [9, 11, 10]],
 [[0, 4, 1], [1, 4, 11], [4, 7, 11], [1, 11, 10]],
 [[10, 1, 4], [1, 3, 4], [3, 8, 4], [11, 10, 4], [7, 11, 4]],
 [[9, 4, 1], [1, 4, 2], [2, 4, 7], [2, 7, 11]],
 [[1, 9, 2], [2, 9, 4], [2, 4, 11], [11, 4, 7], [3, 8, 0]],
 [[11, 4, 7], [2, 4, 11], [0, 4, 2]],
 [[7, 11, 4], [4, 11, 2], [4, 2, 3], [4, 3, 8]],
 [[10, 9, 2], [2, 9, 7], [7, 9, 4], [2, 7, 3]],
 [[2, 10, 7], [10, 9, 7], [9, 4, 7], [0, 2, 7], [8, 0, 7]],
 [[10, 4, 7], [10, 0, 4], [10, 1, 0], [10, 7, 3], [10, 3, 2]],
 [[8, 4, 7], [10, 1, 2]],
 [[4, 1, 9], [7, 1, 4], [3, 1, 7]],
 [[8, 0, 7], [7, 0, 1], [7, 1, 9], [7, 9, 4]],
 [[0, 7, 3], [0, 4, 7]],
 [[8, 4, 7]],
 [[9, 8, 10], [10, 8, 11]],
 [[3, 11, 0], [0, 11, 9], [9, 11, 10]],
 [[0, 10, 1], [8, 10, 0], [11, 10, 8]],
 [[11, 10, 3], [10, 1, 3]],
 [[1, 9, 2], [2, 9, 11], [11, 9, 8]],
 [[9, 2, 1], [11, 2, 9], [3, 11, 9], [0, 3, 9]],
 [[8, 2, 0], [8, 11, 2]],
 [[11, 2, 3]],
 [[2, 8, 3], [10, 8, 2], [9, 8, 10]],
 [[0, 2, 9], [2, 10, 9]],
 [[3, 2, 8], [8, 2, 10], [8, 10, 1], [8, 1, 0]],
 [[1, 2, 10]],
 [[3, 1, 8], [1, 9, 8]],
 [[9, 0, 1]],
 [[3, 0, 8]],
 []]

# Function to perform Marching Cubes on a single cell of the grid.
def marching_cubes_3d_single_cell(f, x, y, z):
    """
    Perform the Marching Cubes algorithm on a single cube cell.

    Parameters:
    - f: Function to evaluate the scalar field at a given point.
    - x, y, z: Coordinates of the lower-left corner of the cube cell.

    Returns:
    - Mesh: A Mesh object containing the vertices and triangles for the isosurface in the cell.
    """
    # Evaluate the function f at each vertex of the cube to determine if the vertex is inside or outside the isosurface.
    f_eval = [None] * 8
    for v in range(8):
        v_pos = VERTICES[v]
        f_eval[v] = f(x + v_pos[0] * CELL_SIZE,
                      y + v_pos[1] * CELL_SIZE,
                      z + v_pos[2] * CELL_SIZE)

    # Determine the index of the case based on which vertices are inside the isosurface.
    case = sum(2**v for v in range(8) if f_eval[v] > 0)

    # Retrieve the faces (triangles) corresponding to this case from the lookup table.
    faces = cases[case]

    def edge_to_boundary_vertex(edge):
        """
        Interpolate to find the vertex at the midpoint of the specified edge.

        Parameters:
        - edge: Index of the edge.

        Returns:
        - V3: The interpolated vertex position as a V3 object.
        """
        v0, v1 = EDGES[edge]
        f0 = f_eval[v0]
        f1 = f_eval[v1]
        t0 = CELL_SIZE - adapt(f0, f1)  # Calculate interpolation factor for vertex 0.
        t1 = CELL_SIZE - t0            # Calculate interpolation factor for vertex 1.
        vert_pos0 = VERTICES[v0]
        vert_pos1 = VERTICES[v1]
        return V3(x + vert_pos0[0] * t0 + vert_pos1[0] * t1,
                  y + vert_pos0[1] * t0 + vert_pos1[1] * t1,
                  z + vert_pos0[2] * t0 + vert_pos1[2] * t1)

    output_verts = []  # List to hold the vertices of the resulting mesh.
    output_tris = []   # List to hold the triangles of the resulting mesh.

    for face in faces:
        # For each face (triangle) in the case, determine the vertices and create triangles.
        edges = face
        verts = list(map(edge_to_boundary_vertex, edges))
        next_vert_index = len(output_verts) + 1
        tri = Tri(
            next_vert_index,
            next_vert_index+1,
            next_vert_index+2,
        )
        output_verts.extend(verts)
        output_tris.append(tri)
    return Mesh(output_verts, output_tris)


def marching_cubes_3d(f, xmin=XMIN, xmax=XMAX, ymin=YMIN, ymax=YMAX, zmin=ZMIN, zmax=ZMAX):
    """
    Applies the Marching Cubes algorithm to a 3D grid within a specified range.

    Parameters:
    - f: Function to evaluate the scalar field.
    - xmin, xmax, ymin, ymax, zmin, zmax: Bounds of the 3D grid.

    Returns:
    - Mesh: A Mesh object representing the isosurface in the specified 3D region.
    """
    mesh = Mesh()  # Initialize an empty Mesh object to store the final result.

    # Iterate over the grid cells in the specified range.
    for x in frange(xmin, xmax, CELL_SIZE):
        for y in frange(ymin, ymax, CELL_SIZE):
            for z in frange(zmin, zmax, CELL_SIZE):
                cell_mesh = marching_cubes_3d_single_cell(f, x, y, z)
                mesh.extend(cell_mesh)
    return mesh


def circle_function(x, y, z):
    """
    Scalar field function representing a sphere.

    Parameters:
    - x, y, z: Coordinates.

    Returns:
    - float: Value of the scalar field at (x, y, z).
    """
    return 2.5 - math.sqrt(x*x + y*y + z*z)


def make_circle_obj(filename):
    """
    Writes an obj file containing a sphere mesh generated by the Marching Cubes algorithm.

    Parameters:
    - filename: Name of the output obj file.
    """
    mesh = marching_cubes_3d(circle_function)
    with open(filename, "w") as f:
        make_obj(f, mesh)


def make_cases_obj():
    """
    Writes obj files demonstrating the main cases of the Marching Cubes algorithm.

    Generates obj files for visualizing various cases to aid in understanding the algorithm.
    """
    import marching_cubes_gen as gen
    assert CELL_SIZE == 1  # Ensure that the cell size is 1 for this demonstration.

    mesh = Mesh()  # Initialize an empty Mesh object to store the cases.
    highlights = Mesh()  # Initialize an empty Mesh object for highlighting vertices.
    offset = V3(0, 0, 0)  # Offset for positioning the cases in the output file.

    # Iterate over all base cases from the Marching Cubes generator.
    for bits, faces in sorted(gen.BASE_CASES.items()):
        verts = set(gen.bits_to_verts(bits))

        # Create a function for this specific case where the scalar field is inside or outside the isosurface.
        def f(x, y, z):
            vert = gen.VERTICES.index((x, y, z))
            return 1 if vert in verts else -1

        # Generate the mesh for this case and add it to the main mesh.
        case_mesh = marching_cubes_3d_single_cell(f, 0, 0, 0, cell_size=1)
        case_mesh = case_mesh.translate(offset)
        mesh.extend(case_mesh)

        # Highlight the solid vertices for visualization.
        highlight = Mesh([V3(*gen.VERTICES[v]) for v in verts], [])
        highlight = highlight.translate(offset)
        highlights.extend(highlight)

        # Update the offset for the next case.
        offset.x += 1.5
        if offset.x > 6:
            offset.x = 0
            offset.y += 1.5

    # Write the generated cases and highlights to obj files.
    with open("cases.obj", "w") as f:
        make_obj(f, mesh)
    with open("case_highlights.obj", "w") as f:
        make_obj(f, highlights)

__all__ = ["marching_cubes_3d"]

if __name__ == "__main__":
    # Generate an obj file for a sphere and save it.
    make_circle_obj("output.obj")
    # Uncomment the line below to generate obj files demonstrating the Marching Cubes cases.
    # make_cases_obj()
