from abc import ABC, abstractmethod

class InvalidRingForMethod(Exception):
    pass

"""Abstract base class for search methods."""

"""
AStar and GreedyBFS called with (heuristic, RING) --> dets for which rings in ring.py
MA_Search called with no args
"""

class SearchClass(ABC):

    @abstractmethod
    def go(self, start):
        pass

    def search(self, start):

        seq, V = self.go(start)
        display("Path found with sequence ",seq)
        display("Remaining cost-free matrix ",V)

        return seq, V

    def check_sequence(self,V,sequence,mat): 

        #gens_dict must be a dictionary of generators with g.symb as key
        #search through gens_dict for key matching start of sequence and post-multiply to recover U
        
        while sequence != "I-":
            for key in self.gens.keys():
                if sequence.startswith(key):
                    sequence = sequence[len(key):]
                    V=V*self.gens[key].mat

        display(V)
        return V == mat