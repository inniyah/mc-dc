This is demonstration code to accompany a <a href="http://www.boristhebrave.com/2018/04/15/marching-cubes-tutorial/">series</a> <a href="http://www.boristhebrave.com/2018/04/15/marching-cubes-3d-tutorial/">of</a> <a href="http://www.boristhebrave.com/2018/04/15/dual-contouring-tutorial/">articles</a> on meshing. Please consult those
articles for more details.

It covers how to implement 2d / 3d Marching Cubes and Dual Contouring. 

There's also additional code to nicely render the results, but that is not polished for re-use.
 
# Usage
 
Simply import one of `marching_cubes_2d.marching_cubes_2d`, `marching_cubes_3d.marching_cubes_3d`,
`dual_contour_2d.dual_contour_2d`, `dual_contour_3d.dual_contour_3d`.

Each function takes an evaluation function, `f`, that determines whether a point is inside or outside
by returning a positive or negative number. The Dual Contouring functions take an additional argument, 
`f_normal`, that returns the gradient as a `V2` or `V3` object. You can optionally pass the range of values
to evaluate f over. The cell size is always 1.

The 2d meshing functions return a unordered list of `common.Edge` objects, the 3d ones return a `utils_3d.Mesh` object.

# License

[CC0]([https://wiki.creativecommons.org/wiki/CC0)

CC0 is the "no copyright reserved" option in the Creative Commons toolkit - it effectively means relinquishing all copyright and similar rights that you hold in a work and dedicating those rights to the public domain.

CC0 is a single purpose tool, designed to take on the dedication function of the former, deprecated Public Domain Dedication and Certification.

How effectively CC0 works will depend on the legal regime in which the work is used, but the tool is intended to effectively release rights even in jurisdictions where it is difficult to do so.

Note that CC0 is a three-tier instrument. We recognize that a waiver may not be effective in some jurisdictions. CC0's enforceability is not solely dependent on the waiver. The fall back public license -- the second tier -- is similar to our Attribution-only license but without the attribution requirement. The third tier is a non-assertion by the copyright holder that even if the waiver and license do not operate as intended, the copyright holder will not take any actions that prevent a user of the work from exercising rights consistent with the intention of the copyright holder as expressed in CC0.

# Introduction to Marching Cubes and Dual Contouring

Marching Cubes and Dual Contouring are algorithms used in computer graphics and computational geometry for surface extraction from scalar fields. These techniques are particularly useful for visualizing 3D medical data, such as MRI scans, or for generating complex 3D models in procedural generation contexts.
1. Marching Cubes (2D and 3D)

Marching Cubes is a classical algorithm used to extract a polygonal mesh of an isosurface from a scalar field (a grid of values). The idea is to iterate over the grid and for each cell (in 2D, a square; in 3D, a cube), determine how the isosurface intersects it.

    2D Marching Cubes:
        In 2D, the algorithm is sometimes called Marching Squares. The grid consists of squares, and each square's corners can either be inside or outside the isosurface.
        The algorithm checks the values at each corner of a square. Based on these values (positive or negative), the square is classified into one of 16 possible configurations (since each corner can be in or out).
        For each configuration, the algorithm determines where the isosurface intersects the square's edges and then connects these points to form line segments.

    3D Marching Cubes:
        Extending the idea to 3D, the grid is made up of cubes, and each cube has eight corners.
        There are 256 possible configurations for how the isosurface can intersect a cube.
        The algorithm determines where the isosurface intersects the edges of each cube and connects these points to form triangles, which together approximate the surface.

Advantages:

    Marching Cubes is relatively simple and easy to implement.
    It generates a high-quality mesh that approximates the surface well.

Disadvantages:

    The generated meshes can be quite large because Marching Cubes does not optimize for mesh simplicity.

2. Dual Contouring (2D and 3D)

Dual Contouring is an advanced technique that, like Marching Cubes, generates a mesh from a scalar field but with several key differences that address some of Marching Cubes' shortcomings.

    Concept:
        Dual Contouring works by placing vertices not at the edges of grid cells but at points within the cells themselves, specifically at the points where the isosurface is likely to pass.
        It tries to minimize the error between the scalar field and the resulting mesh by placing vertices where they best approximate the surface.

    Process:
        For each cell in the grid (square in 2D, cube in 3D), the algorithm calculates where the isosurface is most likely to pass by using the gradients (normals) of the scalar field at the cell's corners.
        A vertex is placed in the cell, and then edges (in 2D) or faces (in 3D) are created by connecting these vertices.

Advantages:

    Dual Contouring tends to produce much simpler meshes than Marching Cubes because it places vertices more intelligently.
    It can produce better results, especially for surfaces with sharp features or corners.

Disadvantages:

    Dual Contouring is more complex to implement because it requires additional information (the gradient of the scalar field).
