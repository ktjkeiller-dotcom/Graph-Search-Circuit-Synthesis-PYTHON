# raised when an invoked matrix operation cannot be completed because the
# dimensions of the two matrices are incompatible. For example, in order to add
# two matrices their dimensions must be identical.
class IncompatibleMatrixSizesException(Exception):
    pass


# raised when you invoke a matrix operation that only be completed on square
# matrices on a non-square matrix.
class MatrixNotSquareException(Exception):
    pass

#raised when trying to pass a fractional matrix into matrix
class EntriesNotIntegerException(Exception):
    pass


# raised when you try and create a matrix with an invalid number of specified
# elements. For example, if you tried to create a 2x2 matrix while only
# specifying 3 elements, this exception would be raised.
class InvalidNumberOfElements(Exception):
    pass


# This class represents an integer matrix over a PID
class Matrix:

    # This hashmap saves solutions to determinant subproblems that have been
    # saved on this matrix to speed up determinant calculations via dynamic
    # programming. For details of its usage see the determinant method.
    _determinant_saves = {}

    # h: the height of the matrix]
    # w: the width of the matrix
    # elements: a list of h*w elements of a PID
    def __init__(self, h, w, content):
        if (isinstance(content,list)):
            if len(content) != h * w:
                raise InvalidNumberOfElements
            self.elements = content

        self.h = h
        self.w = w
        self.elementT = type(content[0])


    # return the sum of the two matrices If the matrices have sizes that
    # prohibit them from being summed a IncompatibleMatrixSizesException is
    # raised
    def __add__(self, y):
        if self.h != y.h or self.w != y.w:
            raise IncompatibleMatrixSizesException

        newElements = []
        for i in range(self.h * self.w):
            newElements.append(self.elements[i] + y.elements[i])
        return type(self)(self.h, self.w, newElements)

    # return the product of the two matrices If the matrices have sizes that
    # prohibit them from being multiplied a IncompatibleMatrixSizesException
    # is raised
    def __mul__(self, y):
        if self.w != y.h:
            raise IncompatibleMatrixSizesException

        newH = self.h
        newW = y.w
        newElements = []
        for i in range(newH):
            for j in range(newW):
                newElement = self.elements[0].get_zero()
                for k in range(self.w):
                    newElement += (self.get(i, k) * y.get(k, j))
                newElements.append(newElement)
        return type(self)(newH, newW, newElements)

    def __matmul__(self, y):
        newH = self.h * y.h
        newW = self.w * y.w
        newElements = []

        for i in range(self.h):
            for bi in range(y.h):
                for j in range(self.w):
                    for bj in range(y.w):
                        newElements.append(self.get(i, j) * y.get(bi, bj))

        return type(self)(newH, newW, newElements)

    # returns a string representation of the matrix
    def __str__(self):
        rows = []
        for i in range(self.h):
            row = self.elements[i * self.w:(i + 1) * self.w]
            rows.append("[ " + " ¦ ".join(map(str, row)) + " ]")
        return "\n".join(rows)

    def _repr_latex_(self):
        rows = []

        for i in range(self.h):
            row = self.elements[i*self.w:(i+1)*self.w]
            rows.append(" & ".join(x.latex() for x in row))

        body = r" \\ ".join(rows)

        return rf"\begin{{bmatrix}}{body}\end{{bmatrix}}"


    # returns whether two matrices are equal, that is whether they have
    # equivalent sizes and are filled with identical elements at each positions
    def __eq__(self, y):
        if self.h != y.h or self.w != y.w:
            return False
        for i in range(self.w * self.h):
            if self.elements[i] != y.elements[i]:
                return False
        return True

    # returns whether two matrices are not equal to one another
    def __ne__(self, y):
        return not self == y

    # returns the determinant of the matrix. This algorithm works on a matrix
    # of size NxN by considering the determinant of some of the submatrices
    # of size (N-1)x(N-1). To speed up this operation we use dynamic
    # programming to memorize the determinant values of submatrices once they
    # have been calculated.

    def determinant(self):
        if self.h != self.w:
            raise MatrixNotSquareException()

        if self.h == 1:
            return self.get(0, 0)

        elements_tuple = tuple(self.elements)
        if elements_tuple in Matrix._determinant_saves:
            return Matrix._determinant_saves[elements_tuple]

        total = type(self.get(0, 0)).get_zero()
        for i in range(self.h):
            scale = self.get(i, 0)
            if i % 2 == 1:
                scale = -scale
            subcontent = []
            for j in range(self.h):
                if i == j:
                    continue
                else:
                    for k in range(1, self.h):
                        subcontent.append(self.get(j, k))
            subMatrix = Matrix(self.h - 1, self.h - 1, subcontent)
            total += scale * subMatrix.determinant()

        Matrix._determinant_saves[elements_tuple] = total
        return total

    # returns the identity matrix of height and width dim, with elements of
    # type elementType
    @staticmethod
    def id(dim, elementType):
        elements = [elementType.get_zero() for _ in range(dim * dim)]
        for i in range(dim):
            elements[i * dim + i] = elementType.get_one()
        return Matrix(dim, dim, elements)

    # returns the matrix element in row i and column j
    def get(self, i, j):
        return self.elements[i * self.w + j]

    # sets the matrix element in row i and column j to be e
    def set(self, i, j, e):
        self.elements[i * self.w + j] = e

    # returns a copy of the matrix
    def copy(self):
        return type(self)(self.h, self.w, self.elements[:])

    def make_integer_mat(self):
        return self.copy()

