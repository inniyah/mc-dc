#!/usr/bin/env python3

"""Provides a function for performing 2D Dual Contouring"""

import math

from common import Edge, adapt, frange
from settings import ADAPTIVE, XMIN, XMAX, YMIN, YMAX, CELL_SIZE
from utils_2d import V2, make_svg
from qef import solve_qef_2d


def dual_contour_2d_find_best_vertex(f, f_normal, x, y):
    """
    Finds the optimal vertex within a grid cell using the Dual Contouring method.

    Parameters:
    - f: A function that computes the scalar field value at a given (x, y) position.
    - f_normal: A function that computes the normal vector of the scalar field at a given (x, y) position.
    - x: The x-coordinate of the lower-left corner of the grid cell.
    - y: The y-coordinate of the lower-left corner of the grid cell.

    Returns:
    - V2: The optimal 2D vertex (position) within the cell, or None if the cell does not contain a sign change.
    """
    if not ADAPTIVE:
        # If not using adaptive refinement, return the center of the cell.
        return V2(x + 0.5 * CELL_SIZE, y + 0.5 * CELL_SIZE)

    # Evaluate the scalar field function at the corners of the cell.
    x0y0 = f(x + 0.0, y + 0.0)
    x0y1 = f(x + 0.0, y + CELL_SIZE)
    x1y0 = f(x + CELL_SIZE, y + 0.0)
    x1y1 = f(x + CELL_SIZE, y + CELL_SIZE)

    # For each edge, identify where there is a sign change
    changes = []
    if (x0y0 > 0) != (x0y1 > 0):
        changes.append([x + 0, y + adapt(x0y0, x0y1)])
    if (x1y0 > 0) != (x1y1 > 0):
        changes.append([x + CELL_SIZE, y + adapt(x1y0, x1y1)])
    if (x0y0 > 0) != (x1y0 > 0):
        changes.append([x + adapt(x0y0, x1y0), y + 0])
    if (x0y1 > 0) != (x1y1 > 0):
        changes.append([x + adapt(x0y1, x1y1), y + CELL_SIZE])

    if len(changes) <= 1:
        # If there are fewer than 2 sign changes, the cell does not contribute to the contour.
        return None

    # For each sign change location v[i], we find the normal n[i].
    normals = []
    for v in changes:
        n = f_normal(v[0], v[1])
        normals.append([n.x, n.y])

    # Solve the quadratic error function to find the best vertex.
    v = solve_qef_2d(x, y, changes, normals)

    return v

def dual_contour_2d(f, f_normal, xmin=XMIN, xmax=XMAX, ymin=YMIN, ymax=YMAX):
    """
    Computes the 2D dual contouring of a scalar field over a specified range of cells.

    Parameters:
    - f: A function that computes the scalar field value at a given (x, y) position.
    - f_normal: A function that computes the normal vector of the scalar field at a given (x, y) position.
    - xmin: The minimum x-coordinate of the range of cells.
    - xmax: The maximum x-coordinate of the range of cells.
    - ymin: The minimum y-coordinate of the range of cells.
    - ymax: The maximum y-coordinate of the range of cells.

    Returns:
    - List of Edge objects representing the boundaries of the contour.
    """

    # Iterates over a cells of size one between the specified range, and evaluates f and f_normal to produce
    # a boundary by Dual Contouring. Returns an unordered list of Edge objects.

    # For each cell, find the the best vertex for fitting f

    # Find the best vertex for each cell.
    verts = {}
    for ix, x in enumerate(frange(xmin, xmax, CELL_SIZE)):
        for iy, y in enumerate(frange(ymin, ymax, CELL_SIZE)):
            verts[ix, iy] = dual_contour_2d_find_best_vertex(f, f_normal, x, y)

    # For each cell edge, emit an edge between the center of the adjacent cells if it is a sign changing edge

    # Create edges where sign changes occur.
    edges = []
    # Do all the vertical sign changes
    for ix, x in enumerate(frange(xmin, xmax, CELL_SIZE)):
        if ix == 0:
            continue
        for iy, y in enumerate(frange(ymin, ymax, CELL_SIZE)):
            y0 = y
            y1 = y + CELL_SIZE
            y0_solid = f(x, y0) > 0
            y1_solid = f(x, y1) > 0
            if y0_solid != y1_solid:
                edges.append(Edge(verts[ix - 1, iy], verts[ix, iy]).swap(y0_solid))
    # Do all the horizontal sign changes
    for iy, y in enumerate(frange(ymin, ymax, CELL_SIZE)):
        if iy == 0:
            continue
        for ix, x in enumerate(frange(xmin, xmax, CELL_SIZE)):
            x0 = x
            x1 = x + CELL_SIZE
            x0_solid = f(x0, y) > 0
            x1_solid = f(x1, y) > 0
            if x0_solid != x1_solid:
                edges.append(Edge(verts[ix, iy - 1], verts[ix, iy]).swap(x0_solid))
    return edges


def circle_function(x, y):
    """
    Defines a circular scalar field centered at the origin with a radius of 2.5.

    Parameters:
    - x: The x-coordinate.
    - y: The y-coordinate.

    Returns:
    - float: The scalar field value at the given coordinates.
    """
    return 2.5 - math.sqrt(x*x + y*y)


def circle_normal(x, y):
    """
    Computes the normal vector to the circular scalar field.

    Parameters:
    - x: The x-coordinate.
    - y: The y-coordinate.

    Returns:
    - V2: The normal vector at the given coordinates.
    """
    l = math.sqrt(x*x + y*y)
    return V2(-x / l, -y / l)


def square_function(x, y):
    """
    Defines a square-shaped scalar field centered at the origin.

    Parameters:
    - x: The x-coordinate.
    - y: The y-coordinate.

    Returns:
    - float: The scalar field value at the given coordinates.
    """
    return 2.5 - max(abs(x), abs(y))


def square_normal(x, y):
    """
    Computes the normal vector to the square-shaped scalar field.

    Parameters:
    - x: The x-coordinate.
    - y: The y-coordinate.

    Returns:
    - V2: The normal vector at the given coordinates.
    """
    if abs(x) > abs(y):
        return V2(-math.copysign(x, 1), 0)
    else:
        return V2(0, -math.copysign(y, 1))

def normal_from_function(f, d=0.01):
    """
    Approximates the gradient of a scalar field function f to compute normals.

    Parameters:
    - f: The scalar field function.
    - d: The scale for finite difference approximation. Smaller values yield more accurate results.

    Returns:
    - Function: A function that computes the gradient (normal) of f at given (x, y) coordinates.
    """

    # Given a sufficiently smooth 2d function, f, returns a function approximating of the gradient of f.
    # d controls the scale, smaller values are a more accurate approximation.

    def norm(x, y):
        return V2(
            (f(x + d, y) - f(x - d, y)) / 2 / d,
            (f(x, y + d) - f(x, y - d)) / 2 / d,
        ).normalize()
    return norm


def t_shape_function(x, y):
    """
    Defines a scalar field for a T-shaped contour.

    Parameters:
    - x: The x-coordinate.
    - y: The y-coordinate.

    Returns:
    - float: The scalar field value at the given coordinates, representing a T-shape.
    """
    if (x, y) in ((0, 0), (0, 1), (0, -1), (1, 0)):
        return 1
    return -1


def intersect_function(x, y):
    """
    Defines a scalar field for an intersection of two shapes.

    Parameters:
    - x: The x-coordinate.
    - y: The y-coordinate.

    Returns:
    - float: The scalar field value at the given coordinates, representing the intersection shape.
    """
    y -= 0.3
    x -= 0.5
    x = abs(x)
    # x += x*x / 1000
    return min(x - y, x + y)


__all__ = ["dual_contour_2d"]

if __name__ == "__main__":
    # Generate contour edges for a circular shape and save as SVG file.
    edges = dual_contour_2d(circle_function, circle_normal)
    with open("example.svg", "w") as file:
        make_svg(file, edges, circle_function)
