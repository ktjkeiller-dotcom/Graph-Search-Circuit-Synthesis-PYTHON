Notes on gate sets!!

CYCLOTOMIC:
- This gate set covers matrices with entries in the fractional ring Z[w, 1/sqrt2]
- We have access to the matrices {I, X, Z, S, H, T, B} 

    1 Qubit:
    - The free set is {S,H} and the cost set is {T}
    - We get 3 generators

    2 Qubits:
    - If we pick C = {CX, H, S} and G = {T}, we get 15 generators however this cannot synthesise matrices with determinant != 1 (without an ancilla). https://arxiv.org/abs/1908.06076 
    - if we pick G = {T, CT} we get 42 generators. This allows to make the matrices with det != 1 and can be used for the GreedyBFS search.

GAUSSIAN:
- This gate set covers matrices with entries in the fractional ring Z[i, 1/2]
- We have access to the matrices: {I, X, Z, S, H, B} (where B is complex, and H is phase adjusted to be rational)

    2 Qubits:
    - We take C = {CX, H, S} and G = {CS}
    - We get 15 generators

ZSQRT:
- This gate set covers matrices where the entries are in the fractional ring Z[1/sqrt2]
- We have access to {I, X, Z, H, B} (where B is irrational and real)

    2 Qubits:
    - We take C = {H, X, Z, CZ} and G = {CH}
    - We get 15 generators

DYADIC: #incomplete
- This gate set covers the 3 qubit matrices where the entries are in the fractional ring Z[1/2]
- We have access to {I, X, Z, HH, B} (where B is for the rational subgroup of the clifford group https://arxiv.org/abs/2404.17677)

    3 Qubits:
    - We take C = {X, CX, HH} and and G = {CCZ, CCX} (verified that these individually are integer in the b-basis)
    - We have to use two different bases to evaluate. In this case, the integer basis is the union of the automorphism groups of two different lattice bases
    - In one basis, CCX is integer, in the other CCZ is integer

