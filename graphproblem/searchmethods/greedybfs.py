from graphproblem.vertex import Vertex
from graphproblem.heuristics.vbheuristic import vb
from graphproblem import searchclass
from graphproblem import generator
from graphproblem.ringenum import Ring

"""This is the search method described in Algorithm 1 of https://doi.org/10.48550/arXiv.2405.19302

This method assumes that at every vertex, there is a generator such that vb vector of
the node its edge leads to strictly weakly majorises that of the current vertex. 

This does not guarantee finding the optimal T-count circuit, however it is guaranteed to terminate
and quickly as it never backtracks.

It relies on the "intermediate lattice property" which was exhaustively demonstrated for the lattices
in the 2-qubit cases of ZI, Z8 (where the costly set is {T, CT}), ZSQRT2 and uses the vb vector defined in vbheuristic.py 
"""

class GreedyBFS(searchclass.SearchClass):
    rings = [Ring.CYCLOTOMIC2_with_CT, Ring.ZSQRT2_2, Ring.GAUSSIAN2] #Rings from Table 4.1 of https://doi.org/10.48550/arXiv.2405.19302

    def __init__(self,RING):

            if RING not in GreedyBFS.rings:
                raise searchclass.InvalidRing

            free_set, cost_set, basis = RING.build_problem_gates()
            self.basis = basis
            self.free = free_set
            self.gens = generator.make_gen_set(None, self.basis, cost_set, free_set)

    def go(self,start):
        #returns string sequence of generators' symbols, and leftover clifford, C
        #start = C * sequence of generators

        mat = self.basis[0]*start*self.basis[1]
        node = Vertex(mat)
        node.sequence="I-"
        node.h = vb(node.mat)
        display(node.mat)

        while (node.mat).denom_exp != 0:

            print("Expanding node with sequence: "+node.sequence)
            display(node.h)

            for generator in self.gens.values():

                #tries different generators until it finds one which strictly weakly majorises
                #then ignores the rest of the generators and move on to the node that strictly weakly majorises

                new_matrix = node.mat * generator.bmat
                neighbour = Vertex(new_matrix)
                neighbour.h = vb(neighbour.mat)
                neighbour.sequence = generator.symb + node.sequence
                print("Created node with sequence: ",neighbour.sequence)
                print(neighbour.h)

                if strictly_majorise(neighbour.h,node.h): #check condition for "strictly weakly majorising"
                    print("Majorised node found")
                    node = neighbour
                    break

        return node.sequence, self.basis[1]*node.mat*self.basis[0]

def strictly_majorise(score1, score2): #returns true if score1 strictly weakly majorises score2

    assert(len(score1)==len(score2))
    sum1=0
    sum2=0
    min_index=False

    for i in range(0,len(score1)):
        sum1+=score1[i]
        sum2+=score2[i]

        if sum1<sum2 and min_index ==False:
            min_index=True
        if sum2<sum1 and min_index == True: 
            min_index=False

    return min_index