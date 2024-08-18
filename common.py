"""Contains utility functions and classes that are shared across different meshing methods."""

import settings

class Edge:
    def __init__(self, v1, v2):
        """
        Initializes an Edge object with two vertices, v1 and v2.

        Parameters:
        v1, v2: These are the vertices that define the endpoints of the edge.
        Vertices are typically represented as tuples or coordinate pairs in a 2D or 3D space.
        """
        self.v1 = v1
        self.v2 = v2

    def swap(self, swap=True):
        """
        Returns a new Edge object with its vertices swapped if the swap parameter is True.

        Parameters:
        swap (bool): If True, returns a new Edge with vertices v1 and v2 swapped.
                     If False, returns a new Edge with the original vertex order.

        Returns:
        Edge: A new Edge object with vertices possibly swapped.
        """
        if swap:
            return Edge(self.v2, self.v1)
        else:
            return Edge(self.v1, self.v2)

def adapt(v0, v1):
    """
    Calculates the interpolation factor between two scalar values, v0 and v1, needed to reach zero.

    This function is used to determine the exact point where the isosurface (represented by the zero level set)
    intersects the edge between the two scalar values, v0 and v1. The function assumes that v0 and v1 have
    opposite signs, meaning the isosurface crosses between them.

    Parameters:
    v0, v1: Scalar values at the endpoints of an edge. They must have opposite signs, meaning one is positive
            and the other is negative.

    Returns:
    float: The distance from v0 to the intersection point (where the scalar field value is zero) along the edge.
           The result is scaled by the CELL_SIZE defined in the settings.
           If ADAPTIVE is False, it returns the midpoint of the edge.
    """
    assert (v1 > 0) != (v0 > 0), "v0 and v1 must have opposite signs for interpolation to make sense."
    if settings.ADAPTIVE:
        return (0 - v0) / (v1 - v0) * settings.CELL_SIZE
    else:
        return 0.5 * settings.CELL_SIZE

def frange(start, stop, step=1):
    """
    Generates a range of floating-point numbers from start to stop, incrementing by step.

    This function is similar to Python's built-in range() but works with floating-point numbers.
    It is useful in situations where you need to iterate over a sequence of non-integer values.

    Parameters:
    start (float): The starting value of the sequence.
    stop (float): The end value of the sequence (the last value generated will be less than stop).
    step (float): The increment between consecutive values in the sequence. Default is 1.

    Yields:
    float: The next value in the sequence from start to stop with the specified step size.
    """
    v = start
    while v < stop:
        yield v
        v += step
