import heapq
from graphproblem.vertex import Vertex
from graphproblem.searchclass import SearchClass
from graphproblem import generator

"""The A* algorithm uses a heuristic function to choose which node in the problem graph to expand next.

 -- The next node to expand has the lowest sum of f = h (heuristic valuation) + g (cost to get to vertex from start) 
  and is guaranteed to find the shortest path when the heuristic is consistent. 
 -- cost to get to vertex from start is the sum of the costs of the generators used 
 -- We can change the weightings in the sum h (heuristic) + g (actual cost) but this does come at the price
  of not guaranteeing a shortest path. Increasing the weighting of the heuristic valuation leads to a more guided
  search and terminates faster, although may not be the optimal circuit. 
  """

class AStar(SearchClass):

    def __init__(self, heuristic, RING):
        
        free_set, cost_set, basis = RING.build_problem_gates()

        self.heuristic = heuristic
        self.basis = basis
        self.free = free_set
        self.gens = generator.make_gen_set(heuristic, self.basis, cost_set, free_set)
    
    def go(self, start): 

        #RETURNS string sequence representing the circuit and left-over clifford, C
        #start = C * generators represented by string sequence

        start_node = Vertex(self.basis[0]*start*self.basis[1]) #start vertex matrix stored in bw basis
        start_node.h = self.heuristic.evaluate(start_node.mat)
        start_node.f = start_node.h
        print("Made start node")

        open_list = []
        closed_set = set()
        open_set = {start}

        heapq.heappush(open_list, start_node)
        #the front of the queue holds the unexpanded node with the lowest f = h + g value and is the most
        #promising node to expand next

        while open_list:

            current_node = heapq.heappop(open_list)
            open_set.discard(current_node.mat)
            closed_set.add(current_node.mat)

            print("Expanding node with path: ", current_node.sequence)
            print("Current node has scores: h(",current_node.h,") + g(",current_node.g,") = f(",current_node.f,")")

            #test that current node isn't clifford. if it is, we are done.
            if current_node.h == 0:
                return current_node.sequence, self.basis[1]*current_node.mat*self.basis[0] #return in normal basis

            for generator in self.gens.values():

                #for each generator, go down that edge to a neighbour
                #calculate f = h + g for neighbour and add it to open_list
                #heapq automatically sorts the unexpanded vertices in open_list according to their f value

                new_matrix = current_node.mat * generator.bmat
                new_node = Vertex(new_matrix)

                if new_node.mat in closed_set:

                    ##currently only checking if the exact matrix is already in the unexpanded node set. 
                    #to properly check if this is a unique vertex, we should use the hnf of the matrix.
                    #***
                    
                    continue

                new_node.h = self.heuristic.evaluate(new_matrix)
                new_node.g = current_node.g + generator.cost
                new_node.f = new_node.g + new_node.h
                
                new_node.sequence = generator.symb + current_node.sequence

                print("New node created with sequence", new_node.sequence)
                print("New node score: h(",new_node.h,") + g(",new_node.g,") = f(",new_node.f,")")
                #display(new_node.mat)

                if new_node.mat not in open_set:

                    #same as above ***

                    heapq.heappush(open_list, new_node)
                    open_set.add(new_node.mat)

        return None
