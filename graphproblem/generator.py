from matrixnormalforms import hnfproblem
from collections import deque
from matrixnormalforms import matrix, qmatrix

"""The generators represent the edges in the problem graph.
-- Each generator's matrix (self.mat) is made from a the matrix of an element from the 'cost_set' 
conjugated by a sequence of elements from the 'free_set'
-- To move through the graph, we post-multiply the vertex' matrix (in the BW basis) by the inverse
of the generator's matrix (in the BW basis), until we end up with just a Clifford matrix.
"""

class Generator:

    def __init__(self, mat,mat_inv, symb,symb_inv):

        self.mat = mat
        self.mat_inv = mat_inv

        self.symb = symb #a string which represents the generator. Eg. H-T-H' is 1 qubit operator representing multiplication H*T*H
        self.symb_inv = symb_inv #represents the inverse of the generator. Eg. H-T-H'
                
        self.bmat = 0 #inverse matrix in the BW basis --> what we post-multiply by to traverse the graph
        self.cost = 0 #the valuation of the generator using the heuristic function
        self.hnf = 0 #used to compare potential generators to ensure we only keep ones with a unique hnf
        

def make_gen_set(heuristic,basis,g_set,free_set): #RETURNS dictionary of generators

    #given the free_set, cost_set, create the finite set of generators with unique Hermite Normal Form 
    #which will represent the edges of the synthesis problem graph

    #breadth first graph search to find the generating set from section 3.2 of https://doi.org/10.48550/arXiv.2405.19302
    #vertex = unique generators
    #edges = conjugation of the matrix representing a unique generator by an element of the free_set. eg. T -> HTH'

    #we expand nodes and discard neighbours that don't represent unique generators (HNFs)

    vertices = [] #stores hnf_problem objects representing each unique vertex/ generator
    working_set = deque() #stores unexpanded nodes

    for g in g_set: 

        #automatically add all the generators in the given cost_set to finite generating set
        #compute their hnf and representations in BW basis

        g.bmat = basis[0]*g.mat*basis[1]
        mat = g.bmat.copy()
        if (isinstance(mat, qmatrix.QMatrix)):
            mat = mat.make_integer()
        assert (isinstance(mat, matrix.Matrix))
        h = hnfproblem.HNFProblem(mat)
        h.compute_hnf()
        g.hnf = h.J
        vertices.append(h)
        working_set.append(g)
        print("Generator added: ",g.symb)
        display(g.mat)

    while working_set: 

        #while there are unexpanded nodes which represent found generators keep looking
        #search ends when every end is a 'dead end' and no more unique HNFs can be found.

        print(len(vertices))
        g = working_set.popleft() #remove node from unexpanded node set

        for c in free_set:

            #to expand a node, test every neighbour by conjugating with each free_set element and 
            #add neighbour to list of unique vertices and to working set, if the HNF is unique.

            if len(vertices)==16:
                print("Manually stopping generator search")
                break

            new_gen = c.conjugate(g)
            new_gen.bmat = basis[0]*new_gen.mat*basis[1]
            mat = new_gen.bmat.copy()
            if (isinstance(mat, qmatrix.QMatrix)):
                mat = mat.make_integer()
            assert (isinstance(mat, matrix.Matrix))
            h = hnfproblem.HNFProblem(mat)
            h.compute_hnf()
            unique = True

            for v in vertices: #compare with every unique generator already found
                if v.hnf_equivalent(h):
                    unique= False
                    print("HNF of ",new_gen.symb," is not unique.")
                    break

            if unique: #we have found a new generator, so add to list of vertices and list of nodes to be expanded
                new_gen.hnf = h.J
                vertices.append(h)
                working_set.append(new_gen)
                g_set.append(new_gen)
                print("Generator added: ",new_gen.symb)
                display(new_gen.mat)

    for g in g_set: 

        #once found the unique generating set, evaluate their costs and bmats (self.inv in the BW basis)

        g.bmat = basis[0]*g.mat_inv*basis[1]

        if heuristic:
            g.cost = heuristic.evaluate(g.bmat)
            
        display(g.symb)
        display(g.mat)

    gen_dict = dict((i.symb,i) for i in g_set) #create a dictionary where the 
    #key is the string symbol representing the generator

    return gen_dict

