"""Contains utilities common to 3D meshing methods, including basic vector operations, 
triangles, quadrilaterals, and mesh handling."""

import math

class V3:
    """A vector in 3D space, representing a point or direction.

    Attributes:
    x (float): The x-coordinate of the vector.
    y (float): The y-coordinate of the vector.
    z (float): The z-coordinate of the vector.
    """
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def normalize(self):
        """Normalize the vector to unit length.

        This method scales the vector so that its magnitude becomes 1, but its direction remains the same.

        Returns:
        V3: A new normalized vector.
        """
        d = math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
        return V3(self.x / d, self.y / d, self.z / d)


class Tri:
    """A 3D triangle defined by three vertices.

    Attributes:
    v1 (V3): The first vertex of the triangle.
    v2 (V3): The second vertex of the triangle.
    v3 (V3): The third vertex of the triangle.
    """
    def __init__(self, v1, v2, v3):
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3

    def map(self, f):
        """Apply a function to each vertex of the triangle.

        This method is useful for transforming the vertices, such as by applying a translation or scaling.

        Parameters:
        f (function): A function that takes a vertex as input and returns a transformed vertex.

        Returns:
        Tri: A new triangle with transformed vertices.
        """
        return Tri(f(self.v1), f(self.v2), f(self.v3))


class Quad:
    """A 3D quadrilateral defined by four vertices.

    Attributes:
    v1 (V3): The first vertex of the quadrilateral.
    v2 (V3): The second vertex of the quadrilateral.
    v3 (V3): The third vertex of the quadrilateral.
    v4 (V3): The fourth vertex of the quadrilateral.
    """
    def __init__(self, v1, v2, v3, v4):
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.v4 = v4

    def map(self, f):
        """Apply a function to each vertex of the quadrilateral.

        This method is useful for transforming the vertices, such as by applying a translation or scaling.

        Parameters:
        f (function): A function that takes a vertex as input and returns a transformed vertex.

        Returns:
        Quad: A new quadrilateral with transformed vertices.
        """
        return Quad(f(self.v1), f(self.v2), f(self.v3), f(self.v4))

    def swap(self, swap=True):
        """Optionally reverses the vertex order of the quadrilateral.

        Reversing the order of the vertices changes the orientation of the quadrilateral, 
        which is important in certain contexts such as defining surface normals.

        Parameters:
        swap (bool): If True, the vertices are reversed.

        Returns:
        Quad: The original or reversed quadrilateral, depending on the `swap` parameter.
        """
        if swap:
            return Quad(self.v4, self.v3, self.v2, self.v1)
        else:
            return Quad(self.v1, self.v2, self.v3, self.v4)


class Mesh:
    """A collection of vertices and faces (triangles or quadrilaterals) that form a 3D mesh.

    Attributes:
    verts (list of V3): The vertices of the mesh.
    faces (list of Tri or Quad): The faces of the mesh, each defined by vertices.
    """
    def __init__(self, verts=None, faces=None):
        self.verts = verts or []
        self.faces = faces or []

    def extend(self, other):
        """Extend the current mesh by adding vertices and faces from another mesh.

        This method is useful for combining multiple meshes into one.

        Parameters:
        other (Mesh): The mesh to add to the current mesh.
        """
        l = len(self.verts)
        f = lambda v: v + l
        self.verts.extend(other.verts)
        self.faces.extend(face.map(f) for face in other.faces)

    def __add__(self, other):
        """Define the addition operator for Mesh objects, allowing for mesh concatenation.

        This method creates a new mesh that is the result of adding two meshes together.

        Parameters:
        other (Mesh): The mesh to add to the current mesh.

        Returns:
        Mesh: A new mesh that combines the vertices and faces of both meshes.
        """
        r = Mesh()
        r.extend(self)
        r.extend(other)
        return r

    def translate(self, offset):
        """Translate (move) the entire mesh by a given offset.

        This method creates a new mesh where each vertex is translated by the given vector.

        Parameters:
        offset (V3): The vector by which to translate the mesh.

        Returns:
        Mesh: A new mesh with translated vertices.
        """
        new_verts = [V3(v.x + offset.x, v.y + offset.y, v.z + offset.z) for v in self.verts]
        return Mesh(new_verts, self.faces)


def make_obj(f, mesh):
    """Export the mesh to a Wavefront OBJ file format.

    This function writes the vertices and faces of the mesh to a file in the OBJ format,
    which is a standard format for 3D models.

    Parameters:
    f (file-like object): The file to write the OBJ data to.
    mesh (Mesh): The mesh to export.
    """
    for v in mesh.verts:
        f.write("v {} {} {}\n".format(v.x, v.y, v.z))
    for face in mesh.faces:
        if isinstance(face, Quad):
            f.write("f {} {} {} {}\n".format(face.v1, face.v2, face.v3, face.v4))
        if isinstance(face, Tri):
            f.write("f {} {} {}\n".format(face.v1, face.v2, face.v3))
