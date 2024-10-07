#!/usr/bin/env python3

import math

# Configuration for the algorithm
ADAPTIVE = True
XMIN = -3
XMAX = 3
YMIN = -3
YMAX = 3
CELL_SIZE = 1

class V2:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def normalize(self):
        d = math.sqrt(self.x*self.x+self.y*self.y)
        return V2(self.x / d, self.y / d)

    def __repr__(self):
        return "{}({})".format(
            type(self).__name__,
            ', '.join([ ("{}='{}'"  if isinstance(v, str) else "{}={}").format(k, v) \
                for k, v in self.__dict__.items() ])
        )

class Edge:
    def __init__(self, v1, v2, prop=None):
        self.v1 = v1  # Start vertex of the edge
        self.v2 = v2  # End vertex of the edge
        self.prop = prop  # Property value associated with the edge

    def swap(self, swap=True):
        if swap:
            return Edge(self.v2, self.v1, self.prop)
        else:
            return Edge(self.v1, self.v2, self.prop)

def adapt(v0, v1):
    """v0 and v1 are numbers of opposite sign. This returns how far you need to interpolate from v0 to v1 to get to 0."""
    assert (v1 > 0) != (v0 > 0), "v0 and v1 do not have opposite sign"
    if ADAPTIVE:
        return (0 - v0) / (v1 - v0) * CELL_SIZE
    else:
        return 0.5 * CELL_SIZE

def frange(start, stop, step=1):
    """Like range, but works for floats"""
    v = start
    while v < stop:
        yield v
        v += step

def make_svg(file, edges, f, f_property=None, palette=None):
    def element(e, **kwargs):
        """Utility function used for rendering SVG elements"""
        s = "<" + e
        for key, value in kwargs.items():
            s += " {}='{}'".format(key, value)
        s += "/>\n"
        return s

    """Writes an SVG file showing the given edges and their properties"""
    scale = 50
    file.write("<?xml version='1.0' encoding='UTF-8'?>\n")
    file.write("<svg version='1.1' xmlns='http://www.w3.org/2000/svg' viewBox='{} {} {} {}'>\n".format(
        XMIN*scale, YMIN*scale, (XMAX-XMIN)*scale, (YMAX-YMIN)*scale))

    file.write("<g transform='scale({})'>\n".format(scale))
    # Draw grid
    for x in frange(XMIN, XMAX, CELL_SIZE):
        for y in frange(YMIN, YMAX, CELL_SIZE):
            file.write(element("rect", x=x, y=y, width=CELL_SIZE, height=CELL_SIZE,
                               style="stroke: grey; stroke-width: 0.02; fill: none"))

    # Draw filled / unfilled circles
    for x in frange(XMIN, XMAX+CELL_SIZE, CELL_SIZE):
        for y in frange(YMIN, YMAX+CELL_SIZE, CELL_SIZE):
            is_solid = f(x, y) > 0
            if is_solid:
                if f_property is not None:
                    c = palette[f_property(x, y)] or "black"
                else:
                    c = "black"
                fill_color = (c if is_solid else "white")
                file.write(element("circle", cx=x, cy=y, r=0.05,
                               style="stroke: black; stroke-width: 0.04; fill: " + fill_color))
            else:
                file.write(element("circle", cx=x, cy=y, r=0.05,
                               style="stroke: gray; stroke-width: 0.02; fill: white"))

    # Draw edges with different colors based on their property
    for edge in edges:
        # Map the property value to a color (e.g., choose a color from a predefined palette)
        color = palette[edge.prop % len(palette)]
        file.write(element("line", x1=edge.v1.x, y1=edge.v1.y, x2=edge.v2.x, y2=edge.v2.y,
                           style='stroke:{};stroke-width:0.04'.format(color)))

    file.write("</g>\n")
    file.write("</svg>\n")

def split_edge(v1, v2, prop1, prop2):
    """
    Splits an edge into two if the properties differ between the endpoints.
    Returns a list of edges.
    """
    if prop1 == prop2:
        print(f"0> {v1}, {v2}, {prop1}")
        return [Edge(v1, v2, prop1)]
    
    # Calculate the midpoint of the edge
    mid_point = V2((v1.x + v2.x) / 2, (v1.y + v2.y) / 2)
    
    # Determine which property to assign to each new edge
    print(f"1> {v1}, {mid_point}, {prop1}")
    print(f"2> {mid_point}, {v2}, {prop2}")
    
    return [
        Edge(v1, mid_point, prop1),  # Assign the first property to the first half
        Edge(mid_point, v2, prop2)   # Assign the second property to the second half
    ]

def marching_cubes_2d_single_cell(f, f_property, x, y):
    """Returns a list of edges that approximate f's boundary for a single cell, with properties."""
    
    # Evaluate the scalar field
    x0y0 = f(x, y)
    x0y1 = f(x, y + CELL_SIZE)
    x1y0 = f(x + CELL_SIZE, y)
    x1y1 = f(x + CELL_SIZE, y + CELL_SIZE)

    # Evaluate the property field
    p_x0y0 = f_property(x, y)
    p_x0y1 = f_property(x, y + CELL_SIZE)
    p_x1y0 = f_property(x + CELL_SIZE, y)
    p_x1y1 = f_property(x + CELL_SIZE, y + CELL_SIZE)

    case = ((1 if x0y0 > 0 else 0) +
            (2 if x0y1 > 0 else 0) +
            (4 if x1y0 > 0 else 0) +
            (8 if x1y1 > 0 else 0))

    if case == 0 or case == 15:
        return []

    edges = []

    if case == 1 or case == 14:
        a1 = adapt(x0y0, x1y0)
        a2 = adapt(x0y0, x0y1)
        v1 = V2(x + a1, y)
        v2 = V2(x, y + a2)
        p1 = p_x0y0 if a1 <= 0.5 else p_x1y0
        p2 = p_x0y0 if a2 <= 0.5 else p_x0y1
        print(f"Case 1/14: Props: {p1} ({a1:.2f}), {p2} ({a2:.2f})")
        edges += split_edge(v1, v2, p1, p2)

    elif case == 2 or case == 13:
        a1 = adapt(x0y0, x0y1)
        a2 = adapt(x0y1, x1y1)
        v1 = V2(x, y + a1)
        v2 = V2(x + a2, y + CELL_SIZE)
        p1 = p_x0y0 if a1 <= 0.5 else p_x0y1
        p2 = p_x0y1 if a2 <= 0.5 else p_x1y1
        print(f"Case 2/13: Props: {p1} ({a1:.2f}), {p2} ({a2:.2f})")
        edges += split_edge(v1, v2, p1, p2)

    elif case == 4 or case == 11:
        a1 = adapt(x1y0, x1y1)
        a2 = adapt(x0y0, x1y0)
        v1 = V2(x + CELL_SIZE, y + a1)
        v2 = V2(x + a2, y)
        p1 = p_x1y0 if a1 <= 0.5 else p_x1y1
        p2 = p_x0y0 if a2 <= 0.5 else p_x1y0
        print(f"Case 4/11: Props: {p1} ({a1:.2f}), {p2} ({a2:.2f})")
        edges += split_edge(v1, v2, p1, p2)

    elif case == 8 or case == 7:
        a1 = adapt(x0y1, x1y1)
        a2 = adapt(x1y0, x1y1)
        v1 = V2(x + a1, y + CELL_SIZE)
        v2 = V2(x + CELL_SIZE, y + a2)
        p1 = p_x0y1 if a1 <= 0.5 else p_x1y1
        p2 = p_x1y0 if a2 <= 0.5 else p_x1y1
        print(f"Case 8/7: Props: {p1} ({a1:.2f}), {p2} ({a2:.2f})")
        edges += split_edge(v1, v2, p1, p2)

    elif case == 3 or case == 12:
        a1 = adapt(x0y0, x1y0)
        a2 = adapt(x0y1, x1y1)
        v1 = V2(x + a1, y)
        v2 = V2(x + a2, y + CELL_SIZE)
        p1 = p_x0y0 if a1 <= 0.5 else p_x1y0
        p2 = p_x0y1 if a2 <= 0.5 else p_x1y1
        print(f"Case 3/12: Props: {p1} ({a1:.2f}), {p2} ({a2:.2f})")
        edges += split_edge(v1, v2, p1, p2)

    elif case == 5 or case == 10:
        a1 = adapt(x0y0, x0y1)
        a2 = adapt(x1y0, x1y1)
        v1 = V2(x, y + a1)
        v2 = V2(x + CELL_SIZE, y + a2)
        p1 = p_x0y0 if a1 <= 0.5 else p_x0y1
        p2 = p_x1y0 if a2 <= 0.5 else p_x1y1
        print(f"Case 5/10: Props: {p1} ({a1:.2f}), {p2} ({a2:.2f})")
        edges += split_edge(v1, v2, p1, p2)

    elif case == 9:
        a1 = adapt(x0y0, x0y1)
        a2 = adapt(x0y0, x0y1)
        a3 = adapt(x0y1, x1y1)
        a4 = adapt(x1y0, x1y1)
        v1 = V2(x + a1, y)
        v2 = V2(x, y + a2)
        v3 = V2(x + a3, y + CELL_SIZE)
        v4 = V2(x + CELL_SIZE, y + a4)
        p1 = p_x0y0 if a1 <= 0.5 else p_x1y0
        p2 = p_x0y0 if a2 <= 0.5 else p_x0y1
        p3 = p_x0y1 if a3 <= 0.5 else p_x1y1
        p4 = p_x1y0 if a4 <= 0.5 else p_x1y1
        print(f"Case 9: Props1: {p1} ({a1:.2f}), {p2} ({a2:.2f})")
        print(f"Case 9: Props2: {p3} ({a3:.2f}), {p4} ({a4:.2f})")
        edges += split_edge(v1, v2, p1, p2)
        edges += split_edge(v3, v4, p3, p4)

    elif case == 6:
        a1 = adapt(x1y0, x1y1)
        a2 = adapt(x0y0, x1y0)
        a3 = adapt(x0y0, x0y1)
        a4 = adapt(x0y1, x1y1)
        v1 = V2(x + CELL_SIZE, y + a1)
        v2 = V2(x + a2, y)
        v3 = V2(x, y + a3)
        v4 = V2(x + a4, y + CELL_SIZE)
        p1 = p_x1y0 if a1 <= 0.5 else p_x1y1
        p2 = p_x0y0 if a2 <= 0.5 else p_x1y0
        p3 = p_x0y0 if a3 <= 0.5 else p_x0y1
        p4 = p_x0y1 if a4 <= 0.5 else p_x1y1
        print(f"Case 6: Props1: {p1} ({a1:.2f}), {p2} ({a2:.2f})")
        print(f"Case 6: Props2: {p3} ({a3:.2f}), {p4} ({a4:.2f})")
        edges += split_edge(v1, v2, p1, p2)
        edges += split_edge(v3, v4, p3, p4)

    return edges

def marching_cubes_2d(f, f_property):
    """Runs marching cubes on the entire grid defined by XMIN, XMAX, YMIN, YMAX, and CELL_SIZE."""
    edges = []
    for x in frange(XMIN, XMAX, CELL_SIZE):
        for y in frange(YMIN, YMAX, CELL_SIZE):
            edges += marching_cubes_2d_single_cell(f, f_property, x, y)
    return edges

def circle_function(x, y):
    return 2.5 - math.sqrt(x*x + y*y)

def example_property_function(x, y):
    """Example property function that assigns a property value based on position."""
    if x < 0 and y < 0:
        return 0  # Property 0
    elif x >= 0 and y < 0:
        return 1  # Property 1
    elif x < 0 and y >= 0:
        return 2  # Property 2
    else:
        return 3  # Property 3

# Example color palette
palette = ["red", "green", "blue", "yellow"]

if __name__ == "__main__":
    edges = marching_cubes_2d(circle_function, example_property_function)
    with open("example.svg", "w") as file:
        make_svg(file, edges, circle_function, example_property_function, palette)
