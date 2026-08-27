from graphproblem.generator import Generator

"""
This defines the cost_free gates
 -- Which are unimodular matrices in the Barnes-Wall basis as defined in https://doi.org/10.48550/arXiv.2405.19302
 -- The class is so-called since the Clifford-gates are all unimodular matrices in this basis.
 -- The class is similar to Generators, however with added functions
 .conjugate() which pre-multiplies by self.mat and post-multiplies by self.mat_inv
 __mul__ which allows for multiplication of a Clifford x Generator to create Matsumoto-Amano generators HT, SHT

 """

class Clifford:
    def __init__(self,mat,inv,name,inv_name):
        self.mat = mat
        self.mat_inv = inv
        self.symb = name #eg. "H-" The dashes help to identify gates that are in series vs parallel.
        self.symb_inv= inv_name # "H'-" where the ' identifies an inverse.

    def conjugate(self, other): #Returns Generator HTH' from H.conjugate(T), where H is Clifford object, and T is Generator object
        new_mat = self.mat*other.mat*self.mat_inv
        new_name = self.symb+other.symb+self.symb_inv
        new_inv = self.mat*other.mat_inv*self.mat_inv
        new_inv_symb = self.symb+other.symb_inv+self.symb_inv
        return Generator(new_mat,new_inv,new_name,new_inv_symb) 

    def __mul__(self, g): #g is a Generator. H*T returns Generator object representing (HT)
        new_mat = self.mat*g.mat
        new_name = self.symb+ g.symb
        new_inv = g.mat_inv * self.mat_inv
        new_inv_symb = g.symb_inv + self.symb_inv
        return Generator(new_mat, new_inv, new_name, new_inv_symb)
