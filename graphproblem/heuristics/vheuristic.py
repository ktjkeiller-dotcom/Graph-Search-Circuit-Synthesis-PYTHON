from graphproblem.heuristicclass import HeuristicClass

class VHeuristic(HeuristicClass):

    #basic heuristic = denominator exponent of DELTA for the ring. 

    def evaluate(self, bmat):
        return bmat.denom_exp