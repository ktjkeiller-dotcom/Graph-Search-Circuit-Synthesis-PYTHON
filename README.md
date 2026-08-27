# Graph-Search-Synthesis-PYTHON

Given a unitary matrix, U, with its entries in a defined Eucliean Domain, algorithm finds its decomposition into a sequence of generators made from "costly" gates and "cost-free" gates.

The generators are made from one "costly" gate and sequences of "cost-free" gates. 

For example, in the Clifford + T gate set, the "costly" gate is T, and the generators are made from "CTC^-1" where C is a sequence of Clifford gates.

Graph traversal methods include two A* searches with different heuristics and Greedy-BFS. The vertices represent matrices that are the same up to pre-multiplication by a "cost-free" matrix and the edges represent generators. Finding the shortest path through the graph corresponds to finding the decomposition which uses the fewest "costly" gates.

Matrix normal forms include Smith Normal Form (for an A* VBheuristic) and Hermite Normal Form (to create the generating sets, and can also be used to show vertices are equivalent/ unique).
