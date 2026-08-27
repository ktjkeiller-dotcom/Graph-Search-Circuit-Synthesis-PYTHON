from matrixnormalforms.hnfproblem import HNFProblem

class Vertex:
    """Represents a vertex in the synthesis graph

    -- Each vertex represents an equivalence class of matrices.

    -- Matrices A,B are equivalent if there is a unimodular matrix C such that A = CB (same HNF)
    Since Cliffords are unimodular in BW basis, this means the matrices have the same cost-optimal sequence
    of generators.
"""

    def __init__(self, mat):
        self.g = 0  # actual cost to get to this vertex from the start vertex
        self.h = 0  # heuristic score of this vertex (measure of how close it is to a solution vertex)
        self.f = 0  # total score (dependent on the type of search implemented. In case of A*, f = g + h

        self.mat = mat #stored in the b-basis
        self.sequence = "I-" #default representation of a vertex is the symbol for the identity matrix
        self.hnf_problem = 0

    def __eq__(self, other): 

        if self.hnf_problem == 0:
            p = HNFProblem(self.mat)
            p.compute_hnf()
            self.hnf_problem = p

        if other.hnf_problem == 0:
            p = HNFProblem(other.mat)
            p.compute_hnf()
            other.hnf_problem = p

        return self.hnf_problem.hnf_equivalent(other.hnf_problem)

    def __lt__(self, other):
        return self.f < other.f  # compare scores. used by heapq to determine what order to expand nodes

    def __hash__(self):
        return hash("com.ktjkeiller.graphproblem",self.f)