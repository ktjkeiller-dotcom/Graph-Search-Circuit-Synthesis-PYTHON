Notation and variable naming:

DELTA refers to the epsilon in the epsilon-rings. I have defined it for each ED according to Table 2.3
of "Multiqubit synthesis and hermitian lattices"

Gates naming conventions:
"H-" refers to the one-qubit Hadamard gate
"HT-IS-" refers to I then H acting on qubit 1, and S then T acting on qubit 2
"CIZ-" is the CZ gate where qubit 1 is the control, and qubit 3 is the target. qubit 2 is left alone.
"S'-" with apostrophe is the inverse of the s gate
The dashes help to signify where a gates are applied in series vs in parallel.

- Vertex matrices are always in the BW basis. 
- Generators have .mat, .mat_inv, .bmat --> .bmat is usually .mat_inv in the BW basis (this is what the A* search in "Multiqubit synthesis and hermitian lattices" needs)
- Cliffords have .mat, .mat_inv which are NOT in the BW basis