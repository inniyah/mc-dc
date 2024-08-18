#!/usr/bin/env python3

"""Provides a function for performing 3D Dual Contouring"""

from common import adapt, frange
from settings import ADAPTIVE, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, CELL_SIZE
import numpy as np
import math
from utils_3d import V3, Tri, Quad, Mesh, make_obj
from qef import solve_qef_3d

# Define this to use triangles instead of quadrilaterals.
USE_TRI = False

def dual_contour_3d_find_best_vertex(f, f_normal, x, y, z):
    """
    Finds the optimal vertex within a 3D grid cell using the Dual Contouring method.

    Parameters:
    - f: A function that computes the scalar field value at a given (x, y, z) position.
    - f_normal: A function that computes the normal vector of the scalar field at a given (x, y, z) position.
    - x: The x-coordinate of the lower-left corner of the grid cell.
    - y: The y-coordinate of the lower-left corner of the grid cell.
    - z: The z-coordinate of the lower-left corner of the grid cell.

    Returns:
    - V3: The optimal 3D vertex (position) within the cell, or None if the cell does not contain a sign change.
    """
    if not ADAPTIVE:
        # If not using adaptive refinement, return the center of the cell.
        return V3(x + 0.5 * CELL_SIZE, y + 0.5 * CELL_SIZE, z + 0.5 * CELL_SIZE)

    # Evaluate the scalar field function at each corner of the cell.
    v = np.empty((2, 2, 2))
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                v[dx, dy, dz] = f(x + dx * CELL_SIZE, y + dy * CELL_SIZE, z + dz * CELL_SIZE)

    # Detect sign changes along the edges of the cell.
    # There are 4 edges along each of the three axes.
    changes = []
    for dx in (0, 1):
        for dy in (0, 1):
            if (v[dx, dy, 0] > 0) != (v[dx, dy, 1] > 0):
                changes.append((x + dx * CELL_SIZE, 
                                y + dy * CELL_SIZE,
                                z + adapt(v[dx, dy, 0], v[dx, dy, 1])))

    for dx in (0, 1):
        for dz in (0, 1):
            if (v[dx, 0, dz] > 0) != (v[dx, 1, dz] > 0):
                changes.append((x + dx * CELL_SIZE,
                                y + adapt(v[dx, 0, dz], v[dx, 1, dz]),
                                z + dz * CELL_SIZE))

    for dy in (0, 1):
        for dz in (0, 1):
            if (v[0, dy, dz] > 0) != (v[1, dy, dz] > 0):
                changes.append((x + adapt(v[0, dy, dz], v[1, dy, dz]),
                                y + dy * CELL_SIZE,
                                z + dz * CELL_SIZE))

    if len(changes) <= 1:
        # If there are fewer than 2 sign changes, the cell does not contribute to the contour.
        return None

    # Compute normals at each sign change location.
    # Solve the quadratic error function to find the best vertex.
    # The error term we are minimizing is || A * x - b ||^2, where A and b are matrices derived from v and n.

    normals = []
    for v in changes:
        n = f_normal(v[0], v[1], v[2])
        normals.append([n.x, n.y, n.z])

    # Solve the quadratic error function in 3D to determine the best vertex position.
    return solve_qef_3d(x, y, z, changes, normals)


def dual_contour_3d(f, f_normal, xmin=XMIN, xmax=XMAX, ymin=YMIN, ymax=YMAX, zmin=ZMIN, zmax=ZMAX):
    """
    Computes the 3D dual contouring of a scalar field over a specified range of cells.

    Parameters:
    - f: A function that computes the scalar field value at a given (x, y, z) position.
    - f_normal: A function that computes the normal vector of the scalar field at a given (x, y, z) position.
    - xmin: The minimum x-coordinate of the range of cells.
    - xmax: The maximum x-coordinate of the range of cells.
    - ymin: The minimum y-coordinate of the range of cells.
    - ymax: The maximum y-coordinate of the range of cells.
    - zmin: The minimum z-coordinate of the range of cells.
    - zmax: The maximum z-coordinate of the range of cells.

    Returns:
    - Mesh: A Mesh object containing the vertices and faces representing the contour.
    """
    # Find the best vertex for each cell.
    vert_array = []
    vert_indices = {}
    for ix, x in enumerate(frange(xmin, xmax, CELL_SIZE)):
        for iy, y in enumerate(frange(ymin, ymax, CELL_SIZE)):
            for iz, z in enumerate(frange(zmin, zmax, CELL_SIZE)):
                vert = dual_contour_3d_find_best_vertex(f, f_normal, x, y, z)
                if vert is None:
                    continue
                vert_array.append(vert)
                vert_indices[ix, iy, iz] = len(vert_array)

    # Create faces where sign changes occur.
    faces = []
    for ix, x in enumerate(frange(xmin, xmax, CELL_SIZE)):
        for iy, y in enumerate(frange(ymin, ymax, CELL_SIZE)):
            for iz, z in enumerate(frange(zmin, zmax, CELL_SIZE)):
                if x > xmin and y > ymin:
                    solid1 = f(x, y, z + 0) > 0
                    solid2 = f(x, y, z + CELL_SIZE) > 0
                    if solid1 != solid2:
                      if not USE_TRI:
                        faces.append(Quad(
                            vert_indices[(ix - 1, iy - 1, iz)],
                            vert_indices[(ix - 0, iy - 1, iz)],
                            vert_indices[(ix - 0, iy - 0, iz)],
                            vert_indices[(ix - 1, iy - 0, iz)],
                        ).swap(solid2))
                      else:
                        faces.append(Tri(
                            vert_indices[(ix - 1, iy - 1, iz)],
                            vert_indices[(ix - 0, iy - 1, iz)],
                            vert_indices[(ix - 1, iy - 0, iz)],
                        ).swap(solid2))
                        faces.append(Tri(
                            vert_indices[(ix - 0, iy - 1, iz)],
                            vert_indices[(ix - 0, iy - 0, iz)],
                            vert_indices[(ix - 1, iy - 0, iz)],
                        ).swap(solid2))
                if x > xmin and z > zmin:
                    solid1 = f(x, y + 0, z) > 0
                    solid2 = f(x, y + CELL_SIZE, z) > 0
                    if solid1 != solid2:
                      if not USE_TRI:
                        faces.append(Quad(
                            vert_indices[(ix - 1, iy, iz - 1)],
                            vert_indices[(ix - 0, iy, iz - 1)],
                            vert_indices[(ix - 0, iy, iz - 0)],
                            vert_indices[(ix - 1, iy, iz - 0)],
                        ).swap(solid1))
                      else:
                        faces.append(Tri(
                            vert_indices[(ix - 1, iy, iz - 1)],
                            vert_indices[(ix - 0, iy, iz - 1)],
                            vert_indices[(ix - 1, iy, iz - 0)],
                        ).swap(solid1))
                        faces.append(Tri(
                            vert_indices[(ix - 0, iy, iz - 1)],
                            vert_indices[(ix - 0, iy, iz - 0)],
                            vert_indices[(ix - 1, iy, iz - 0)],
                        ).swap(solid1))
                if y > ymin and z > zmin:
                    solid1 = f(x + 0, y, z) > 0
                    solid2 = f(x + CELL_SIZE, y, z) > 0
                    if solid1 != solid2:
                      if not USE_TRI:
                        faces.append(Quad(
                            vert_indices[(ix, iy - 1, iz - 1)],
                            vert_indices[(ix, iy - 0, iz - 1)],
                            vert_indices[(ix, iy - 0, iz - 0)],
                            vert_indices[(ix, iy - 1, iz - 0)],
                        ).swap(solid2))
                      else:
                        faces.append(Tri(
                            vert_indices[(ix, iy - 1, iz - 1)],
                            vert_indices[(ix, iy - 0, iz - 1)],
                            vert_indices[(ix, iy - 1, iz - 0)],
                        ).swap(solid2))
                        faces.append(Tri(
                            vert_indices[(ix, iy - 0, iz - 1)],
                            vert_indices[(ix, iy - 0, iz - 0)],
                            vert_indices[(ix, iy - 1, iz - 0)],
                        ).swap(solid2))
    # Return a Mesh object containing the vertices and faces.
    return Mesh(vert_array, faces)


def circle_function(x, y, z):
    """
    Defines a spherical scalar field centered at the origin with a radius of 2.5.

    Parameters:
    - x: The x-coordinate.
    - y: The y-coordinate.
    - z: The z-coordinate.

    Returns:
    - float: The scalar field value at the given coordinates.
    """
    return 2.5 - math.sqrt(x*x + y*y + z*z)


def circle_normal(x, y, z):
    """
    Computes the normal vector to the spherical scalar field.

    Parameters:
    - x: The x-coordinate.
    - y: The y-coordinate.
    - z: The z-coordinate.

    Returns:
    - V3: The normal vector at the given coordinates.
    """
    l = math.sqrt(x*x + y*y + z*z)
    return V3(-x / l, -y / l, -z / l)

def intersect_function(x, y, z):
    """
    Defines a scalar field for an intersection of two shapes in 3D.

    Parameters:
    - x: The x-coordinate.
    - y: The y-coordinate.
    - z: The z-coordinate.

    Returns:
    - float: The scalar field value at the given coordinates, representing the intersection shape.
    """
    y -= 0.3
    x -= 0.5
    x = abs(x)
    return min(x - y, x + y)

def normal_from_function(f, d=0.01):
    """
    Approximates the gradient of a scalar field function f to compute normals.

    Parameters:
    - f: The scalar field function.
    - d: The scale for finite difference approximation. Smaller values yield more accurate results.

    Returns:
    - Function: A function that computes the gradient (normal) of f at given (x, y, z) coordinates.
    """
    def norm(x, y, z):
        return V3(
            (f(x + d, y, z) - f(x - d, y, z)) / 2 / d,
            (f(x, y + d, z) - f(x, y - d, z)) / 2 / d,
            (f(x, y, z + d) - f(x, y, z - d)) / 2 / d,
        ).normalize()
    return norm

__all__ = ["dual_contour_3d"]

if __name__ == "__main__":
    # Generate a 3D mesh for a spherical shape and save it as an OBJ file.
    mesh = dual_contour_3d(circle_function, normal_from_function(circle_function))
    with open("output.obj", "w") as f:
        make_obj(f, mesh)
