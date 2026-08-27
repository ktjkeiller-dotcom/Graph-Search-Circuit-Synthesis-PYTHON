from matrixnormalforms import matrix, snfproblem, q, qmatrix, ed
from graphproblem.heuristicclass import HeuristicClass

class VBHeuristic(HeuristicClass): 

    #Heuristic defined in Equation 25 with pk = delta(k,n)

    def evaluate(self, bmat):
        vb_vector = vb(bmat)
        total = 0

        for i in vb_vector:
            total += i

        return total

def vb(mat):

    #calculate the vb(U) vector from Equation 22

    #smith normal form is only defined for integer matrices but mat is not necessarily an integer matrix 
    #so we multiply by the largest denominator to make it integer
    #this does not affect its normal form being unique

    assert(mat.h==mat.w and mat.h%2==0) #matrix must be square and have even dimension
    if (isinstance(mat, qmatrix.QMatrix)):
        exp = mat.denom_exp
        delta =q.Q(mat.int_type.DELTA,0) #need to make fractional version of integer DELTA so we can multiply
        for i in range(0,exp):
            mat = mat*delta 
        mat = mat.make_integer()
    else:
        assert (isinstance(mat, matrix.Matrix))
        exp=0

    prob = snfproblem.SNFProblem(mat)
    prob.compute_snf()

    diagonal = [prob.J.elements[i*prob.J.w +i] for i in range(0,int(prob.J.h/2))] #First half of diagonal entries 

    assert((isinstance(v,ed.Ed) for v in diagonal)) #elements in diagonal should be integers
    delta = prob.J.elementT.DELTA #the integer version of DELTA
    vb_vector = []

    for entry in diagonal:
        k=0 #power
        while entry.is_divisible(delta):
            k+=1
            z = delta.matrix.LUsolve(entry.coeffs)
            entry = prob.J.elementT([int(i) for i in z])

        vb_vector.append(exp-k) #effectively undoing the multiplications we did earlier to make it an integer matrix

    vb_vector.reverse() 
    #the diagonal elements of snf form divide the next but we need vb_vector entries to be non-increasing so
    #we permute the elements of the diagonal. this is allowed since permutations are unimodular.

    return vb_vector