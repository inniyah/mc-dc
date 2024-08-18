"""
Global settings for meshing algorithms (Marching Cubes and Dual Contouring).
These settings help control various aspects of the algorithms, reducing the need
for repetitive argument passing. They also allow you to explore different variants
of the algorithms, as discussed in the accompanying articles.
"""

# ADAPTIVE controls whether the algorithms adaptively select vertices based on
# the underlying function's characteristics. If set to False, the midpoint of each
# cell is used instead, simplifying the mesh generation process. This is mainly
# for illustration purposes to show the difference between adaptive and non-adaptive meshing.
ADAPTIVE = True

# In Dual Contouring, the CLIP setting forces the selected vertex to remain within
# the cell boundaries. When set to True, this can lead to less accurate but more
# contained mesh vertices, which might be useful for certain applications.
CLIP = False

# The BOUNDARY setting in Dual Contouring applies constraints during the vertex
# minimization process, ensuring the vertex does not stray outside the cell during
# optimization. This is crucial for maintaining the integrity of the mesh when sharp
# features or boundaries are present.
BOUNDARY = True

# BIAS adds additional penalties during the vertex optimization process in Dual Contouring.
# This encourages the algorithm to keep the vertex within the cell, which can prevent
# vertices from drifting too far from their expected positions. It is particularly useful
# when the mesh needs to conform closely to the original grid.
BIAS = True

# BIAS_STRENGTH determines how strong the bias is relative to the gradients used in the
# optimization process. A value of 1.0 would mean the bias is as strong as the input gradients.
# Smaller values, like 0.01, make the bias influence weaker, allowing the mesh to still
# adapt but with some preference for staying within the cell.
BIAS_STRENGTH = 0.01

# Default bounds for the 3D space over which the scalar field is evaluated. These values
# define the region of interest in the coordinate system. Adjusting these values will
# affect the area of the scalar field that is processed by the meshing algorithms.
XMIN = -3
XMAX = 3
YMIN = -3
YMAX = 3
ZMIN = -3
ZMAX = 3

# The size of each cell in the grid. This defines the resolution of the grid, where
# smaller cell sizes result in finer meshes but require more computational resources.
# Larger cell sizes produce coarser meshes but are less computationally demanding.
CELL_SIZE = 1

# EPSILON is a small value used to prevent floating point precision issues. It is
# particularly important in comparisons or calculations where exact equality of
# floating point numbers could lead to errors due to the nature of floating point arithmetic.
EPS = 1e-8