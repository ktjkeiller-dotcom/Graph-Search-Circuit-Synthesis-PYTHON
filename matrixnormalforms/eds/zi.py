from matrixnormalforms import ed, pid, q
import sympy as sp

"""
This represents the Gaussian integers. 
 -- Every integer in the ring can be written as z = {a + bi : a, b are integers}
 -- Coeffs correspond to [a,b]
 -- Norm is taken to be a^2 + b^2
 """

class ZI(ed.Ed):
    dim = 2 #we require two coefficients to define an integer in this ring

    def __init__(self, content):
        if isinstance(content, list):
            if len(content) != 2:
                raise pid.InvalidInitialContent
            if not (isinstance(content[0], int) and
                    isinstance(content[1], int)):
                raise pid.InvalidInitialContent
            else:
                self.a = content[0]
                self.b = content[1]

        elif isinstance(content, q.Q):
            if content.denom_exp != 0:
                raise pid.InvalidInitialContent
            else:
                self.a = content.integer.a
                self.b = content.integer.b

        elif isinstance(content, ZI):
            self.a = content.a
            self.b = content.b
        else:
            raise pid.InvalidInitialContent

        self.coeffs = sp.Matrix([self.a, self.b])
        self.matrix = sp.Matrix([
            [self.a, -self.b],
            [self.b, self.a]])

    def __str__(self):

        if self.a == 0 and self.b == 0:
            return "0"
        elif self.a == 0:
            return f"{self.b}i"
        elif self.b == 0:
            return f"{self.a}"
        else:
            return f"{self.a}{self.b:+}i"
        
    def _repr_latex_(self):
        basis = ["", r"i"]
        parts = []

        for coeff, name in zip(self.coeffs, basis):
            if coeff == 0:
                continue
            if name == "":
                parts.append(str(coeff))
            elif coeff == 1:
                parts.append(name)
            elif coeff == -1:
                parts.append(f"-{name}")
            else:
                parts.append(f"{coeff}{name}")

        if not parts:
            return "0"
        
        return " + ".join(parts).replace("+ -", "- ")

    def latex(self):
        basis = ["", r"i"]
        parts = []

        for coeff, name in zip(self.coeffs, basis):
            if coeff == 0:
                continue

            if name == "":
                parts.append(str(coeff))
            elif coeff == 1:
                parts.append(name)
            elif coeff == -1:
                parts.append(f"-{name}")
            else:
                parts.append(f"{coeff}{name}")

        if not parts:
            return "0"

        return " + ".join(parts).replace("+ -", "- ")
    
    def __hash__(self):
        return hash(("com.ktjkeiller.matrixnormalforms.zi", self.a, self.b))

    def __eq__(self, y):
        return self.a == y.a and self.b == y.b

    def __ne__(self, y):
        return not self == y

    def __neg__(self):
        return ZI([-self.a, -self.b])

    def __add__(self, y):
        return ZI([self.a + y.a, self.b + y.b])

    def __sub__(self, y):
        return ZI([self.a - y.a, self.b - y.b])

    def __mul__(self, y):
        return ZI([self.a * y.a - self.b * y.b, self.a * y.b + self.b * y.a])

    # return the complex conjugate
    def com(self):
        return ZI([self.a, -self.b])

    # return the product of this ZI with the complex conjugate of another ZI
    def num(self, y):
        return self * y.com()

    #self = r + q * y

    def get_q(self, y):
        n1 = self.num(y).a
        n2 = self.num(y).b
        d = y.a * y.a + y.b * y.b
        comp1 = (n1 + d // 2) // d
        comp2 = (n2 + d // 2) // d
        return ZI([comp1, comp2])

    def get_r(self, y):
        return ZI([(self - y * (self // y)).a, (self - y * (self // y)).b])

    #Returns the square of the standard complex norm

    def norm(self):
        return self.a * self.a + self.b * self.b

    #return the additive identity of the ring
    @staticmethod
    def get_zero():
        return ZI([0, 0])

    # return the multiplicative identity of the ring
    @staticmethod
    def get_one():
        return ZI([1, 0])

    @staticmethod
    def get_scalar(num):
        return ZI([num, 0])


ZI.DELTA = ZI([1,1]) #the denominator for this ring as used in the graph-problem synthesis algorithm

