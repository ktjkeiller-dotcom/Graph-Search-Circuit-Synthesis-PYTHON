from matrixnormalforms import ed

class InvalidInitialContent(Exception):
    pass

class InvalidMultiplicationType(Exception):
    pass


""" This defines a general fractional ring where the numerator belongs to a ring defined as a subclass of ED.
 -- The denominator can take the value of powers of a specific integer, which is defined for each ring as .DELTA
 -- we define basic arithmetic and a simplification check.
 """

class Q():
    def __init__(self, content, de ,int_type = 0): 
        #de = denominator exponent
        #int_type defines the integer ring we are using to create this fractional ring

        if not (isinstance(de,int)): #de must take an integer value
            raise InvalidInitialContent
        
        if isinstance(content, list):
            assert(issubclass(int_type,ed.Ed)) #checks that the int_type has been defined (not default = 0)
            #this allows to define the numerator of the fraction as an element of the ring 'int_type

            if len(content) != int_type.dim: #check we have provided enough coefficients for the numerator ring
                raise InvalidInitialContent
            else:
                self.int_type = int_type
                self.integer = self.int_type(content)

        elif isinstance(content, ed.Ed):
            self.int_type = type(content)
            self.integer = content 

        else:
            raise InvalidInitialContent
        
        #self.integer is the numerator 
        self.denom_exp = de 
        self.denom = self.int_type.DELTA
        self.simplify()

    def __eq__(self, y):
        return (self.integer == y.integer and self.denom == y.denom and self.denom_exp==y.denom_exp)

    def __ne__(self, y):
        return not self==y

    def simplify(self): #reduce self.denom_expo(nent) if self.denom divides self.integer
        while ((self.integer.is_divisible(self.denom)) and (not self.denom_exp ==0)):
            z = self.denom.matrix.LUsolve(self.integer.coeffs)
            z = [int(x) for x in z]
            self.integer = self.int_type(z)
            self.denom_exp -= 1
        return self

    def latex(self):
            if self.denom_exp==0:
                return self.integer.latex()
            else:
                return rf"\frac{{{self.integer.latex()}}}{{{(self.int_type.DELTA.latex())}^{{{self.denom_exp}}}}}"

    def __str__(self):
            """
            Pretty-print the number in LaTeX form.
            """
            if self.denom_exp ==0:
                return rf"${self.integer.latex()}$"
            else:
                return rf"$\frac{{{self.integer.latex()}}}{{{(self.int_type.DELTA.latex())}^{{{self.denom_exp}}}}}$"
            
    def __neg__(self):
        return Q(-self.integer, self.denom_exp)

    def __add__(self, y):
        s = self.copy()
        z = y.copy()
        if (s.denom_exp > z.denom_exp):
            denom = s.denom_exp
            for i in range(0,(s.denom_exp-z.denom_exp)):
                assert(self.denom_exp != 0)
                z.integer = z.integer * self.denom

        elif (s.denom_exp < z.denom_exp):
            denom = z.denom_exp
            for i in range(0,(z.denom_exp-s.denom_exp)):
                s.integer = s.integer * self.denom
        else:
            denom = z.denom_exp
        res = Q(s.integer+z.integer,denom)
        return res.simplify()

    def __sub__(self, y):
        return self + (-y)

    def __mul__(self, y):
        copy = self.copy()
        z=y.copy()
        if (isinstance(z,Q)):
            copy.integer = copy.integer*z.integer
            copy.denom_exp = copy.denom_exp + z.denom_exp
        elif (isinstance(z,self.int_type)):
            copy.integer = copy.integer *z
        else:
            raise InvalidMultiplicationType
        return copy.simplify()

    def copy(self):
        return Q(self.integer.copy(), self.denom_exp)

    def get_zero(self):
        return Q(self.int_type.get_zero(),0)

    def get_one(self):
        return Q(self.int_type.get_one(),0)

    def __hash__(self):
            return hash(("com.ktjkeiller.matrixnormalforms.q", self.integer.dim))