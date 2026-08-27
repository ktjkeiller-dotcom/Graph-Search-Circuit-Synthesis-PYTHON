Matsumoto-amano vs CGC' decompositions

I ran on five tests and both forms gave the same T-count, although looked very different!

I wonder if a two-qubit version of MA_form might work, using the denom_exponent to guide us? In
the case that the CX are not needed. Is there a basis where H,S are integer but CX is not? We could
maybe use that to split the algorithm.

if there are no cx matrices, then the whole matrix is a tensor product of two matsumoto-amano forms, 
and the one-qubit search could be adapted easily.