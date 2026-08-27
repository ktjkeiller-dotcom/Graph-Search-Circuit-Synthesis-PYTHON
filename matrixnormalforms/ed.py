from abc import abstractmethod, ABC
from matrixnormalforms import pid
import sympy as sp

"""
This abstract class represents Euclidean domains
-- Computing the Hermite and Smith normal forms requires the extended gcd algorithm which we define here
   using the get_q and get_r defined uniquely for each ring.
-- We also define .is_unit, .is_divisible, .is_equivalent
"""

class Ed(pid.Pid, ABC):

    def coeffs(self): #integers are stored as coeffients multiplying the basis elements in the ring
        pass 

    def matrix(self): #defines how two integers in the ED multiply. y * x = z --> y.matrix * x.coeffs = z.coeffs
        pass

    @abstractmethod
    def norm(self):
            pass

    # returns whether |self| < |x| (a norm-wise comparison)
    def __lt__(self, x):
        return self.norm() < x.norm()

    # returns whether |self| > |x| (a norm-wise comparison)
    def __gt__(self, x):
        return self.norm() > x.norm()

    # in addition to the norm it is required x.__floordiv__(y) and x.__mod__(y)
    # return q and r respectively such x = q*y + r is a euclidean relation

    @abstractmethod
    def get_q(self, x):
        pass

    def __floordiv__(self, x):
        return self.get_q(x)

    @abstractmethod
    def get_r(self, x):
        pass

    def __mod__(self, x):
        return self.get_r(x)

    # in a euclidean domain, the extended euclidean algorithm can be used
    # to find the extended_gcd, which makes this problem much easier

    def extended_gcd(a, b):
        x0 = type(b).get_one() #multiplicative identity
        x1 = type(b).get_zero() #additive identity
        y0 = type(b).get_zero()
        y1 = type(b).get_one()
        while b != type(b).get_zero(): 

            #while second isn't identically zero
            #standard euclidean algorithm part 

            tempa = a 
            tempb = b

            q = tempa // tempb #floor: calls temps.get_q(tempb) 
            a = tempb 
            b = tempa % tempb #remainder: calls tempa.get_r(tempb)

            #extended part
            tempx0 = x0
            x0 = x1
            x1 = tempx0 - q * x0
            tempy0 = y0
            y0 = y1
            y1 = tempy0 - q * y0

        return [a, x0, y0]

    def is_unit(self):
        if self.norm() == 1:
            return True
        else:
            return False
    
    def is_divisible(self,y): #true if y¦self
        z = y.matrix.LUsolve(self.coeffs)
        return all(isinstance(x, sp.Integer) for x in list(z))
    
    def is_equivalent(self,y): #check if two elements are equivalent up to multiplication by a unit

        res = self.matrix.LUsolve(y.coeffs)
        x = [int(i) for i in res]

        return (type(self)(x).is_unit())
    
    def copy(self):
        return (type(self)(self))
