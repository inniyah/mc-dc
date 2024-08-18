import numpy
import numpy.linalg

from utils_2d import V2
from utils_3d import V3
import settings


class QEF:
    """Represents and solves the Quadratic Error Function (QEF).

    The QEF is used to find the best-fit point for a set of constraints in a least-squares sense.
    It minimizes the error in a quadratic fashion, given by the equation || A * x - b ||^2.
    """

    def __init__(self, A, b, fixed_values):
        """
        Initializes a QEF instance.

        Parameters:
        - A: The matrix representing the coefficients of the linear equations.
        - b: The vector representing the right-hand side of the equations.
        - fixed_values: List of fixed values for certain variables (if any).
        """
        self.A = A
        self.b = b
        self.fixed_values = fixed_values

    def evaluate(self, x):
        """
        Evaluates the QEF at a given point.

        Parameters:
        - x: The point (vector) at which to evaluate the QEF.

        Returns:
        - float: The Euclidean norm of the difference between A*x and b, which is the error to minimize.
        """
        x = numpy.array(x)
        return numpy.linalg.norm(numpy.matmul(self.A, x) - self.b)

    def eval_with_pos(self, x):
        """
        Evaluates the QEF at a position and returns the result in the same format as solve().

        Parameters:
        - x: The position (vector) at which to evaluate the QEF.

        Returns:
        - tuple: A tuple containing the error and the position.
        """
        return self.evaluate(x), x

    @staticmethod
    def make_2d(positions, normals):
        """
        Creates a QEF for a 2D problem.

        Parameters:
        - positions: List of 2D positions (vectors).
        - normals: List of 2D normals (vectors) associated with the positions.

        Returns:
        - QEF: A QEF instance representing the 2D problem.
        """
        A = numpy.array(normals)
        b = [v[0] * n[0] + v[1] * n[1] for v, n in zip(positions, normals)]
        fixed_values = [None] * A.shape[1]
        return QEF(A, b, fixed_values)

    @staticmethod
    def make_3d(positions, normals):
        """
        Creates a QEF for a 3D problem.

        Parameters:
        - positions: List of 3D positions (vectors).
        - normals: List of 3D normals (vectors) associated with the positions.

        Returns:
        - QEF: A QEF instance representing the 3D problem.
        """
        A = numpy.array(normals)
        b = [v[0] * n[0] + v[1] * n[1] + v[2] * n[2] for v, n in zip(positions, normals)]
        fixed_values = [None] * A.shape[1]
        return QEF(A, b, fixed_values)

    def fix_axis(self, axis, value):
        """
        Constrains a specific axis to a fixed value and returns a new QEF.

        Parameters:
        - axis: The index of the axis to be fixed (0 for x, 1 for y, 2 for z).
        - value: The fixed value for the given axis.

        Returns:
        - QEF: A new QEF instance with the specified axis fixed.
        """
        # Adjust the right-hand side vector b to account for the fixed axis.
        b = self.b[:] - self.A[:, axis] * value
        # Remove the fixed axis from matrix A.
        A = numpy.delete(self.A, axis, 1)
        fixed_values = self.fixed_values[:]
        fixed_values[axis] = value
        return QEF(A, b, fixed_values)

    def solve(self):
        """
        Solves the QEF to find the point that minimizes the error.

        Returns:
        - tuple: A tuple containing the residual error and the optimal point (vector).
        """
        result, residual, rank, s = numpy.linalg.lstsq(self.A, self.b, rcond=None)
        if len(residual) == 0:
            residual = self.evaluate(result)
        else:
            residual = residual[0]
        # Construct the final position vector by including fixed values.
        # Result only contains the solution for the unfixed axis,
        # we need to add back all the ones we previously fixed.
        position = []
        i = 0
        for value in self.fixed_values:
            if value is None:
                position.append(result[i])
                i += 1
            else:
                position.append(value)
        return residual, position


def solve_qef_2d(x, y, positions, normals):
    """
    Solves a 2D QEF problem to find the optimal point within a cell.

    Parameters:
    - x, y: Coordinates of the top-left corner of the cell.
    - positions: List of 2D positions (vectors).
    - normals: List of 2D normals (vectors) associated with the positions.

    Returns:
    - V2: The optimal 2D point (vector) within the cell.
    """

    # The error term we are trying to minimize is sum( dot(x-v[i], n[i]) ^ 2)
    # This should be minimized over the unit square with top left point (x, y)

    # In other words, minimize || A * x - b || ^2 where A and b are a matrix and vector
    # derived from v and n
    # The heavy lifting is done by the QEF class, but this function includes some important
    # tricks to cope with edge cases

    # This is demonstration code and isn't optimized, there are many good C++ implementations
    # out there if you need speed.

    CELL_SIZE = settings.CELL_SIZE

    if settings.BIAS:
        # Add extra normals that add extra error the further we go
        # from the cell, this encourages the final result to be
        # inside the cell
        # These normals are shorter than the input normals
        # as that makes the bias weaker,  we want them to only
        # really be important when the input is ambiguous

        # Take a simple average of positions as the point we will
        # pull towards.

        # Add bias normals to encourage the result to stay within the cell.
        mass_point = numpy.mean(positions, axis=0)

        normals.append([settings.BIAS_STRENGTH, 0])
        positions.append(mass_point)
        normals.append([0, settings.BIAS_STRENGTH])
        positions.append(mass_point)

    qef = QEF.make_2d(positions, normals)
    residual, v = qef.solve()

    if settings.BOUNDARY:
        def inside(r):
            return x <= r[1][0] <= x + CELL_SIZE and y <= r[1][1] <= y + CELL_SIZE

        # Check if the solution is within the cell. If not, constrain the QEF to the boundaries.
        if not inside((residual, v)):
            # We constrain the the QEF to the horizontal and vertical
            # lines bordering the cell, and find the best point of those
            r1 = qef.fix_axis(0, x + 0).solve()
            r2 = qef.fix_axis(0, x + CELL_SIZE).solve()
            r3 = qef.fix_axis(1, y + 0).solve()
            r4 = qef.fix_axis(1, y + CELL_SIZE).solve()

            rs = list(filter(inside, [r1, r2, r3, r4]))

            if len(rs) == 0:
                # It's still possible that those lines (which are infinite)
                # cause solutions outside the box. So finally, we evaluate which corner
                # of the cell looks best
                r1 = qef.eval_with_pos((x + 0, y + 0))
                r2 = qef.eval_with_pos((x + 0, y + CELL_SIZE))
                r3 = qef.eval_with_pos((x + CELL_SIZE, y + 0))
                r4 = qef.eval_with_pos((x + CELL_SIZE, y + CELL_SIZE))

                rs = list(filter(inside, [r1, r2, r3, r4]))

            # Pick the best of the available options
            residual, v = min(rs)

    if settings.CLIP:
        # Ensure that the point remains within the cell boundaries.
        v[0] = numpy.clip(v[0], x, x + CELL_SIZE)
        v[1] = numpy.clip(v[1], y, y + CELL_SIZE)

    return V2(v[0], v[1])


def solve_qef_3d(x, y, z, positions, normals):
    """
    Solves a 3D QEF problem to find the optimal point within a cell.

    Parameters:
    - x, y, z: Coordinates of the corner of the cell.
    - positions: List of 3D positions (vectors).
    - normals: List of 3D normals (vectors) associated with the positions.

    Returns:
    - V3: The optimal 3D point (vector) within the cell.
    """

    # The error term we are trying to minimize is sum( dot(x-v[i], n[i]) ^ 2)
    # This should be minimized over the unit square with top left point (x, y)

    # In other words, minimize || A * x - b || ^2 where A and b are a matrix and vector
    # derived from v and n
    # The heavy lifting is done by the QEF class, but this function includes some important
    # tricks to cope with edge cases

    # This is demonstration code and isn't optimized, there are many good C++ implementations
    # out there if you need speed.

    CELL_SIZE = settings.CELL_SIZE

    if settings.BIAS:
        # Add extra normals that add extra error the further we go
        # from the cell, this encourages the final result to be
        # inside the cell
        # These normals are shorter than the input normals
        # as that makes the bias weaker,  we want them to only
        # really be important when the input is ambiguous

        # Take a simple average of positions as the point we will
        # pull towards.

        # Add bias normals to encourage the result to stay within the cell.
        mass_point = numpy.mean(positions, axis=0)

        normals.append([settings.BIAS_STRENGTH, 0, 0])
        positions.append(mass_point)
        normals.append([0, settings.BIAS_STRENGTH, 0])
        positions.append(mass_point)
        normals.append([0, 0, settings.BIAS_STRENGTH])
        positions.append(mass_point)

    qef = QEF.make_3d(positions, normals)
    residual, v = qef.solve()

    if settings.BOUNDARY:
        def inside(r):
            return (x <= r[1][0] <= x + CELL_SIZE and
                    y <= r[1][1] <= y + CELL_SIZE and
                    z <= r[1][2] <= z + CELL_SIZE)

        # Check if the solution is within the cell. If not, constrain the QEF to the cell boundaries.
        if not inside((residual, v)):
            # We constrain the QEF to the 6 planes bordering the cell, and find the best point of those.
            r1 = qef.fix_axis(0, x + 0).solve()
            r2 = qef.fix_axis(0, x + CELL_SIZE).solve()
            r3 = qef.fix_axis(1, y + 0).solve()
            r4 = qef.fix_axis(1, y + CELL_SIZE).solve()
            r5 = qef.fix_axis(2, z + 0).solve()
            r6 = qef.fix_axis(2, z + CELL_SIZE).solve()

            rs = list(filter(inside, [r1, r2, r3, r4, r5, r6]))

            if len(rs) == 0:
                # It's still possible that those planes (which are infinite) cause solutions outside the box.
                # So now try the 12 lines bordering the cell.
                r1  = qef.fix_axis(1, y + 0).fix_axis(0, x + 0).solve()
                r2  = qef.fix_axis(1, y + CELL_SIZE).fix_axis(0, x + 0).solve()
                r3  = qef.fix_axis(1, y + 0).fix_axis(0, x + CELL_SIZE).solve()
                r4  = qef.fix_axis(1, y + CELL_SIZE).fix_axis(0, x + CELL_SIZE).solve()
                r5  = qef.fix_axis(2, z + 0).fix_axis(0, x + 0).solve()
                r6  = qef.fix_axis(2, z + CELL_SIZE).fix_axis(0, x + 0).solve()
                r7  = qef.fix_axis(2, z + 0).fix_axis(0, x + CELL_SIZE).solve()
                r8  = qef.fix_axis(2, z + CELL_SIZE).fix_axis(0, x + CELL_SIZE).solve()
                r9  = qef.fix_axis(2, z + 0).fix_axis(1, y + 0).solve()
                r10 = qef.fix_axis(2, z + CELL_SIZE).fix_axis(1, y + 0).solve()
                r11 = qef.fix_axis(2, z + 0).fix_axis(1, y + CELL_SIZE).solve()
                r12 = qef.fix_axis(2, z + CELL_SIZE).fix_axis(1, y + CELL_SIZE).solve()

                rs = list(filter(inside, [r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12]))

            if len(rs) == 0:
                # Evaluate the corners of the cell if no suitable solution is found on the planes or lines.
                r1 = qef.eval_with_pos((x + 0, y + 0, z + 0))
                r2 = qef.eval_with_pos((x + 0, y + 0, z + CELL_SIZE))
                r3 = qef.eval_with_pos((x + 0, y + CELL_SIZE, z + 0))
                r4 = qef.eval_with_pos((x + 0, y + CELL_SIZE, z + CELL_SIZE))
                r5 = qef.eval_with_pos((x + CELL_SIZE, y + 0, z + 0))
                r6 = qef.eval_with_pos((x + CELL_SIZE, y + 0, z + CELL_SIZE))
                r7 = qef.eval_with_pos((x + CELL_SIZE, y + CELL_SIZE, z + 0))
                r8 = qef.eval_with_pos((x + CELL_SIZE, y + CELL_SIZE, z + CELL_SIZE))

                rs = list(filter(inside, [r1, r2, r3, r4, r5, r6, r7, r8]))

            # Pick the best solution from the available options.
            residual, v = min(rs)

    if settings.CLIP:
        # Ensure that the point remains within the cell boundaries.
        v[0] = numpy.clip(v[0], x, x + CELL_SIZE)
        v[1] = numpy.clip(v[1], y, y + CELL_SIZE)
        v[2] = numpy.clip(v[2], z, z + CELL_SIZE)

    return V3(v[0], v[1], v[2])
