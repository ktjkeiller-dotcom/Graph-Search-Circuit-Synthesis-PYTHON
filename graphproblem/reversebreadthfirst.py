import heapq
from graphproblem import vertex
from graphproblem.heuristics.vbheuristic import vb

"This expands from the identity matrix backwards"
"-- A breadth-first search expands all available nodes before moving on to expand all the newly available nodes"
"-- This makes it a good way to check for loops or patterns in sequences of generators according to the"
"VBHeuristic"""

class rbf_search():

    def __init__(self, gens, basis):
        self.gens = gens
        self.basis = basis

def go(self, start, num_levels):
    #num_levels is maximum number of nodes to expand in any direction
    #it will evaluate (num_generators)^(num_levels) vertices

    start_mat = self.basis[0]*start*self.basis[1] #could use START = identity
    start_node = vertex.Vertex(start_mat)
    start_node.h = vb(start)

    current_set = []
    next_set =[]
    final_set =[]

    for g in self.gens.values(): #redefine g.bmat since this is usually the inverse in the bw basis 
        self.basis[0] * g.mat * self.basis[1] 

    heapq.heappush(current_set,start_node)

    for i in range(0,num_levels): 

        while current_set: #for every node at current level, expand it.

            current_node = heapq.heappop(current_set) #remove this node from current level
            heapq.heappush(final_set,current_node) #add this node to final set of evaluated nodes 
            
            for g in self.gens.values():

                new_matrix = current_node.mat * g.bmat
                new_node = vertex.Vertex(new_matrix)
                new_node.sequence = current_node.sequence + g.symb

                print("New node created with sequence",new_node.sequence)
                heapq.heappush(next_set, new_node) #add all new nodes to next level to be expanded 
                
        current_set = next_set #redefine current_set to next level of nodes 
        next_set = []

    d = dict() 

    while final_set: #sort expanded nodes by their vb scores

        x = heapq.heappop(final_set)
        print("Node with sequence", x.sequence)
        x.h = vb(x.mat)

        if tuple(x.h) in d:
            d[tuple(x.h)].append(x.sequence)
        else:
            d[tuple(x.h)] = [x.sequence]

    while current_set: #sort unexpanded nodes by their vb scores into same buckets as final_set

        y = heapq.heappop(current_set)
        print("Node with sequence", y.sequence)
        y.h = vb(y.mat)

        if tuple(y.h) in d:
            d[tuple(y.h)].append(y.sequence)
        else:
            d[tuple(y.h)] = [y.sequence]

    for key, value in d.items(): #print sorted lists of nodes acording to vb scores 

        print("Score: ",key)
        print("Sequences: ",value)
    
    return d