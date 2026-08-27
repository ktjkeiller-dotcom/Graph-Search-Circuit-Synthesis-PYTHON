from matrixnormalforms import matrix, q

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

"""This defines a matrix where the entries belong to fractional rings as defined in q.py 
 -- additional functionality of this subclass is the attribute .denom_exp which is the 
  greatest denominator exponent of all of the entries in the matrix
"""

class QMatrix(matrix.Matrix):

    def __init__(self, h, w, content):

        #check all elements are from fractional ring
        if (isinstance(content,list) and all(isinstance(content[i],q.Q) for i in range(0,len(content)))):
            self.elements = content
            if (len(content) != h*w):
                raise InvalidNumberOfElements

        if (isinstance(content, QMatrix)): #make a copy
            print(type(content.elements[0]))
            self.elements = content.elements

        elif (isinstance(content, matrix.Matrix)): #turn an integer matrix into a fractional with denom_exp = 0
            print(type(content.elements[0]))
            self.elements = [q.Q(element,0) for element in content.elements]

        self.h=h
        self.w=w
        self.int_type = self.elements[0].int_type
        self.denom_exp = max(element.denom_exp for element in self.elements) 
        #greatest denominator exponent of all entries in the matrix

    def simplify(self): #simplify the fractions and recalculate denom_exp
        self.elements = [element.simplify() for element in self.elements]
        self.denom_exp = max(element.denom_exp for element in self.elements)

    def make_integer(self): ##multiplies every element up by int_type.DELTA^denom_exp to make it an integer matrix
        copy = self.elements.copy()
        for i in range(0,self.denom_exp):
            copy = [element*self.int_type.DELTA for element in copy]
        copy  = [self.int_type(element.integer) for element in copy]
        m =matrix.Matrix(self.h, self.w, copy)
        return m

    def __mul__(self, y): #defines standard matrix multiplication
            if (isinstance(y,q.Q)):
                els = [self.elements[i]*y for i in range(0,len(self.elements))]
                fm = QMatrix(self.h, self.w, els)
                fm.simplify()
                return fm
            else:
                if self.w != y.h:
                    raise matrix.IncompatibleMatrixSizesException
                newH = self.h
                newW = y.w
                newElements = []
                for i in range(newH):
                    for j in range(newW):
                        newElement = self.elements[0].get_zero()
                        for k in range(self.w):
                            newElement += (self.get(i, k) * y.get(k, j))
                        newElements.append(newElement)
                return QMatrix(newH, newW, newElements)

    def __matmul__(self, y): #defines the tensor product
        newH = self.h * y.h
        newW = self.w * y.w
        newElements = []

        for i in range(self.h):
            for bi in range(y.h):
                for j in range(self.w):
                    for bj in range(y.w):
                        newElements.append(self.get(i, j) * y.get(bi, bj))

        return QMatrix(newH, newW, newElements)

    def __hash__(self):
        return hash(("com.ktjkeiller.matrixnormalforms.qmatrix", self.h))
    