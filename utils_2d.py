"""Contains utility functions and classes specific to 2D meshing methods."""

from settings import XMIN, XMAX, YMIN, YMAX, CELL_SIZE, EPS
from common import frange

import math

from typing import List, Callable, Dict, Any, TextIO

class V2:
    def __init__(self, x: float, y: float) -> None:
        """
        Initializes a 2D vector (V2) with x and y coordinates.

        Parameters:
        x (float): The x-coordinate of the vector.
        y (float): The y-coordinate of the vector.
        """
        self.x: float = x
        self.y: float = y

    def normalize(self) -> 'V2':
        """
        Returns a normalized version of the vector (unit vector) that maintains the
        direction of the original vector but scales its length to 1.

        Returns:
        V2: A new V2 object representing the normalized vector.
        """
        d: float = math.sqrt(self.x * self.x + self.y * self.y)  # Calculate the magnitude (length) of the vector
        return V2(self.x / d, self.y / d)  # Divide each component by the magnitude to normalize

def element(e: str, **kwargs: Dict[str, Any]) -> str:
    """
    Utility function to generate an SVG element string with attributes.

    Parameters:
    e (str): The name of the SVG element (e.g., 'rect', 'circle', 'line').
    kwargs (dict): Key-value pairs representing the attributes of the SVG element (e.g., 'x', 'y', 'width', 'height').

    Returns:
    str: A string representing the SVG element with its attributes.
    """
    s: str = "<" + e
    for key, value in kwargs.items():
        s += " {}='{}'".format(key, value)  # Add each attribute to the element string
    s += "/>\n"
    return s

def make_svg(file: TextIO, edges: List['Edge'], f: Callable[[float, float], float]) -> None:
    """
    Writes an SVG file that visually represents the mesh and scalar field.

    This function generates an SVG image showing the grid, the filled and unfilled circles based on the scalar field,
    and the edges of the mesh. The SVG is written directly to the specified file.

    Parameters:
    file (file-like object): The file where the SVG content will be written.
    edges (list of Edge): A list of Edge objects representing the edges of the mesh to be drawn.
    f (function): A scalar field function that takes x, y coordinates and returns a value indicating
                  whether the point is inside or outside the surface (positive or negative value).

    The function works as follows:
    - Draws a grid over the specified coordinate range.
    - Fills circles based on the scalar field value at grid points (solid or unfilled).
    - Draws the edges that form the mesh.
    - Highlights the vertices of the edges with red squares.
    """
    scale: int = 50  # Scale factor for the SVG image to make it visually clear
    file.write("<?xml version='1.0' encoding='UTF-8'?>\n")  # Start of the SVG file with XML declaration
    file.write("<svg version='1.1' xmlns='http://www.w3.org/2000/svg' viewBox='{} {} {} {}'>\n".format(
        XMIN * scale, YMIN * scale, (XMAX - XMIN) * scale, (YMAX - YMIN) * scale))  # Define the SVG canvas

    file.write("<g transform='scale({})'>\n".format(scale))  # Group element to scale the entire SVG content

    # Draw the grid lines within the specified bounds using rectangles
    for x in frange(XMIN, XMAX, CELL_SIZE):
        for y in frange(YMIN, YMAX, CELL_SIZE):
            file.write(element("rect", x=x, y=y, width=CELL_SIZE, height=CELL_SIZE,
                               style="stroke: grey; stroke-width: 0.02; fill: none"))

    # Draw circles at grid points, filled or unfilled based on the scalar field function f
    for x in frange(XMIN, XMAX + CELL_SIZE, CELL_SIZE):
        for y in frange(YMIN, YMAX + CELL_SIZE, CELL_SIZE):
            is_solid: bool = f(x, y) > 0  # Determine if the point is inside the surface
            fill_color: str = "black" if is_solid else "white"  # Black for inside, white for outside
            file.write(element("circle", cx=x, cy=y, r=0.05,
                               style="stroke: black; stroke-width: 0.02; fill: " + fill_color))

    # Draw edges of the mesh as red lines connecting vertices
    for edge in edges:
        file.write(element("line", x1=edge.v1.x, y1=edge.v1.y, x2=edge.v2.x, y2=edge.v2.y,
                           style='stroke:rgb(255,0,0);stroke-width:0.04'))

    # Highlight the vertices of the edges with red squares for better visibility
    r: float = 0.05  # Half of the square's side length
    for v in [v for edge in edges for v in (edge.v1, edge.v2)]:
        file.write(element("rect", x=(v.x - r), y=(v.y - r), width=2 * r, height=2 * r,
                           style='fill: red'))

    file.write("</g>\n")  # Close the group element
    file.write("</svg>\n")  # Close the SVG element
