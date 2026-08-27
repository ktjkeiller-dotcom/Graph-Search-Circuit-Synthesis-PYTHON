from matrixnormalforms import ed, pid, q
import sympy as sp

"""
This represents the standard integer ring. 
-- It extends ED as it is a Euclidean domain, and we provide a remainder and quotient function 
   to be used in the extended Euclidean algorithm.
-- We also define basic operations, a norm function, and the integer which I call .DELTA and is 
   the 'epsilon' if this ring is taken as an 'epsilon-ring'.

https://en.wikipedia.org/wiki/Ring_of_integers
"""

class Z(ed.Ed):
    dim = 1

    def __init__(self, content):
        if isinstance(content, int):
            self.a = content
        elif isinstance(content, list):
            if len(content) != 1 or not isinstance(content[0], int):
                raise pid.InvalidInitialContent
            else:
                self.a = content[0]
        elif isinstance(content, Z):
            self.a= content.a

        elif isinstance(content, q.Q):
            if content.denom_exp != 0:
                raise pid.InvalidInitialContent
            else:
                self.a = content.integer.a

        else:
            raise pid.InvalidInitialContent

        self.coeffs = sp.Matrix([self.a])
        self.matrix = sp.Matrix([[self.a]])

    def _repr_latex_(self):
        basis = [""]
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
            basis = [""]
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
    
    def __str__(self):
        return str(self.a)

    def __hash__(self):
        return hash(("com.ktjkeiller.matrixnormalforms.z", self.a))

    def __eq__(self, x):
        return self.a == x.a

    def __ne__(self, x):
        return self.a != x.a

    def __neg__(self):
        return Z(-self.a)

    def __add__(self, x):
        return Z(self.a + x.a)

    def __sub__(self, x):
        return Z(self.a - x.a)

    def __mul__(self, x):
        return Z(self.a * x.a)


    #self = q + x * r
    def get_q(self, x): 

        return Z(self.a // x.a)

    def get_r(self, x):
        return Z(self.a % x.a)

    # we choose the norm to be the absolute value squared
    def norm(self):
        return self.a * self.a

    # return the additive identity of the ring
    @staticmethod
    def get_zero():
        return Z(0)

    # return the multiplicative identity of the ring
    @staticmethod
    def get_one():
        return Z(1)

    # Returns whether the integer here is a unit. see:
    # https://en.wikipedia.org/wiki/Unit_(ring_theory)
    def is_unit(self):
        return (self.a == 1) or (self.a == -1)

    @staticmethod
    def get_scalar(num):
        return Z([num])

Z.DELTA = Z([2])