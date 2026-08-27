from abc import ABC, abstractmethod

"""A heuristic function estimates how close a vertex is to being a solution vertex
-- We use them in the A* search to guide the search

-- vb_heuristic uses the Smith Normal Form of the matrices 
-- v_heuristic uses the .denom_exp of the vertex' matrix in the b-basis

Details of both heuristics can be found in:
https://doi.org/10.48550/arXiv.2405.19302 "Multi-qubit synthesis and hermitian lattices" - Kliuchnikov, Schonnenbeck
"""

class HeuristicClass(ABC):
    """Abstract base class for heuristics."""

    #not technically needed in python, but i like declaring ABCs

    @abstractmethod
    def evaluate(self, bmat):
        pass
    