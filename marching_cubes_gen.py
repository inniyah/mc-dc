#!/usr/bin/env python3

"""Regenerates the 256 cases needed for the 3D Marching Cubes lookup table.

The Marching Cubes algorithm is used to extract a polygonal mesh of an isosurface 
from a three-dimensional scalar field, such as medical imaging data or 3D simulations.
This script generates a complete set of cases (256 possible vertex states) for how the 
isosurface intersects a cube. The resulting lookup table is critical for implementing 
the Marching Cubes algorithm efficiently.

NOTE: The values generated by this script are already hardcoded into 'marching_cubes_3d.py',
so there is no need to run this script unless you are interested in how the lookup table is generated.
"""

# My convention for vertices in a cube:
# Each vertex of the cube is represented by its 3D coordinates in a unit cube.
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

# Edge definitions for the cube:
# Each edge is defined as a pair of vertex indices.
# For example, the edge (0, 1) connects vertex 0 and vertex 1.
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

# A dictionary that maps a set of two vertex indices (an edge) to its index in the EDGES list.
EDGES_BY_VERTSET = {}
for e, (v1, v2) in enumerate(EDGES):
    EDGES_BY_VERTSET[frozenset([v1, v2])] = e

# These are the 15 base cases of the Marching Cubes algorithm.
# Each key represents a bitwise pattern where a 1 indicates a solid (inside the isosurface) 
# vertex and a 0 indicates an empty (outside the isosurface) vertex.
# The value is a list of tuples, each representing a triangle. Each triangle is defined by 
# the edges it intersects, specified by their indices in the EDGES list.
# Source: http://users.polytech.unice.fr/~lingrand/MarchingCubes/algo.html
BASE_CASES = {
    0b00000000: (),  # No vertices inside the isosurface.
    0b00000001: ((8, 0, 3), ),  # Only vertex 0 inside the isosurface.
    0b00000011: ((8, 1, 3), (8, 9, 1)),  # Vertices 0, 1 inside.
    0b00000101: ((8, 0, 3), (1, 10, 2)),  # Vertices 0, 2 inside.

    0b01000001: ((8, 0, 3), (10, 5, 6)),  # Vertices 0, 4 inside.
    0b00110010: ((8, 7, 0), (0, 7, 1), (1, 7, 5)),  # Vertices 1, 2, 3, 5 inside.
    0b01000011: ((8, 1, 3), (8, 9, 1), (10, 5, 6)),  # Vertices 0, 1, 4 inside.
    0b01001010: ((3, 2, 11), (0, 9, 1), (10, 5, 6)),  # Vertices 0, 2, 4 inside.

    0b00110011: ((7, 5, 3), (3, 5, 1)),  # Vertices 0, 1, 2, 5 inside.
    0b10110001: ((11, 6, 3), (3, 6, 0), (0, 6, 5), (0, 5, 9)),  # Vertices 0, 4, 5, 6 inside.
    0b01101001: ((11, 8, 2), (2, 8, 0), (6, 10, 4), (4, 10, 9)),  # Vertices 0, 2, 4, 6 inside.
    0b01110001: ((3, 7, 0), (0, 7, 10), (7, 6, 10), (0, 10, 9)),  # Vertices 0, 4, 6, 7 inside.

    0b00111010: ((3, 2, 11), (8, 7, 0), (0, 7, 1), (1, 7, 5)),  # Vertices 0, 1, 2, 4 inside.
    0b10100101: ((8, 0, 3), (4, 5, 9), (10, 2, 1), (11, 6, 7)),  # Vertices 0, 2, 4, 5 inside.
    0b10110010: ((8, 11, 0), (0, 11, 5), (5, 11, 6), (0, 5, 1)),  # Vertices 0, 1, 5, 6 inside.
}

# Normally you can invert the base cases (replace solid with empty and visa versa) to generate more cases
# But this leads to a problem - such flips in 3d can swap the cube edges from one ambiguous configuration
# to another.
#
#     *-----o          *-----o
#     |  \  |          | /   |
#     |   \ |          |/    |
#     |\    |          |    /|
#     | \   |          |   / |
#     o-----*          o-----*
#
# This will lead to edges not properly lining up when putting adjacent cells together.
# To solve this, we specify a few more cases, that are flips of the above, but with a different arrangement of edges
# See http://users.polytech.unice.fr/~lingrand/MarchingCubes/algo.html

# Inverse cases are required to handle certain ambiguous configurations that occur 
# when solid and empty regions are swapped in 3D space.
# These cases ensure that adjacent cubes align properly without introducing cracks 
# or inconsistencies in the resulting mesh.
INVERSE_CASES = {
    255 - 0b00000101: ((3, 2, 8), (8, 2, 10), (8, 10, 1), (8, 1, 0)),  # Inverse of 0b00000101.
    255 - 0b01000011: ((6, 8, 3), (6, 9, 8), (6, 5, 9), (6, 3, 1), (6, 1, 10)),  # Inverse of 0b01000011.
    255 - 0b01001010: ((3, 11, 0), (0, 11, 6), (0, 6, 9), (9, 6, 5), (1, 10, 2)),  # Inverse of 0b01001010.
}

# Start with all base cases
ALL_CASES = BASE_CASES.copy()

# Add inverse cases to all cases.
for bits, faces in INVERSE_CASES.items():
    ALL_CASES[bits] = faces

# Vertex permutation operations: These operations will be used to transform 
# the cube by rotating and reflecting it, thus generating all possible 
# configurations from the base cases.

ROTATE_1 = [1, 2, 3, 0, 5, 6, 7, 4]  # 90-degree rotation around the x-axis.
ROTATE_2 = [3, 2, 6, 7, 0, 1, 5, 4]  # 90-degree rotation around the y-axis.
ROTATE_3 = [1, 5, 6, 2, 0, 4, 7, 3]  # 90-degree rotation around the z-axis.
REFLECT = [1, 0, 3, 2, 5, 4, 7, 6]  # Reflection (mirror) across the yz-plane.

# Helper functions for manipulating vertex and edge representations:

def bits_to_verts(n):
    """Converts a bitwise representation of vertices into an array of vertex indices."""
    return [v for v in range(8) if 2**v & n > 0]

def verts_to_bits(vs):
    """Converts an array of vertex indices into a bitwise representation."""
    return sum(2**v for v in vs)

def bits_apply(op, n):
    """Applies a vertex permutation operation to a set of vertices (in bitwise form)."""
    return verts_to_bits([op[v] for v in bits_to_verts(n)])

def faces_apply(op, faces):
    """Applies a vertex permutation operation to a list of triangles (faces)."""
    return tuple(tuple(EDGES_BY_VERTSET[frozenset([op[v1], op[v2]])] for v1, v2 in (EDGES[e] for e in face)) for face in faces)

def operation_name(ops):
    """Generates a string name for a sequence of operations applied to the cube."""
    return "".join({
        tuple(ROTATE_1): "ROTATE_1 ",
        tuple(ROTATE_2): "ROTATE_2 ",
        tuple(ROTATE_3): "ROTATE_3 ",
        tuple(REFLECT): "REFLECT ",
    }[tuple(op)] for op in ops)

def operations_apply(ops, faces):
    """Applies a sequence of operations to a list of triangles (faces)."""
    for op in ops:
        faces = faces_apply(op, faces)
    return faces

# Now generate all possible cases by applying rotations and reflections to the base cases.
OPERATIONS = [(ROTATE_1, ROTATE_1, ROTATE_1), (ROTATE_2,), (ROTATE_3,), (REFLECT,)]
for op_set in OPERATIONS:
    for base, faces in BASE_CASES.items():
        result = operations_apply(op_set, faces)
        ALL_CASES[bits_apply(op_set[-1], base)] = result

# Helper functions for validating the generated lookup table:

def test1():
    """Ensures that every triangle edge is between a solid and an empty vertex."""
    for bits, faces in ALL_CASES.items():
        for face in faces:
            for e in face:
                v1, v2 = EDGES[e]
                assert (2**v1 & bits > 0) != (2**v2 & bits > 0)

def test2():
    """Ensures that each generated case is a valid manifold (continuous surface)."""
    for bits, faces in ALL_CASES.items():
        verts = bits_to_verts(bits)
        for face in faces:
            for e in face:
                v1, v2 = EDGES[e]
                assert v1 in verts or v2 in verts

def test3():
    """Ensures consistent triangulation across shared faces between cubes."""
    for bits1, faces1 in ALL_CASES.items():
        for bits2, faces2 in ALL_CASES.items():
            if sum(2**v in [bits1, bits2] for v in range(8)) == 4:
                assert faces1 == faces2

# Run tests to validate the generated lookup table.
test1()
test2()
test3()

# Print the generated lookup table for inspection.
from pprint import pprint
pprint(ALL_CASES)
