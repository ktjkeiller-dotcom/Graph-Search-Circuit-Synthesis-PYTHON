To use each search method (A*, GreedyBFS, Matsumoto-Amano), instantiate a search_method.
- s = astar.AStar(heuristic, RING) --> A* has two choices of heuristic, and works over different rings
- s = greedybfs.GreedyBFS(RING) --> heuristic/comparison of vertices is specific to Greedy, but can choose ring
- s = ma_search.MASearch() -->no choice, only defined for 1 qubit Clifford+T

Example:

from graphproblem.searchmethods import astar 
from graphproblem.ring import Ring 
from graphproblem.heuristics.vheuristic import VHeuristic

RING = Ring.ZSQRT2_2 --> synthesise a 2-qubit matrix over Z[sqrt2] 
heuristic = LDEHeuristic() --> choose the heuristic to use (LDEHeuristic or KHeuristic) 
search_method = astar.AStar(heuristic, RING) --> choose the search method (A*, GreedyBFS Matsumoto-Amano)

gates = RING.ring_gates 
CH = gates["CH-"] 
EYE = gates["I-"] 
HI = gates["H-"] @ EYE 
IH = EYE @ gates["H-"] 
ZI = gates["Z-"] @ EYE 
CX = gates["CX-"] 
TEST = CH*CX*HI*SI*CH*IH*CX*CH*IH*CH*CX

seq, V = search_method.search(TEST) --> returns string sequence of generators and remaining cost-free gate