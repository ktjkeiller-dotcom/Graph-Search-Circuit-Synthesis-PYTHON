from graphproblem.vertex import Vertex
from graphproblem.searchclass import SearchClass
from graphproblem.generator import Generator
from graphproblem.ringenum import Ring

"""
Matsumoto-amano form decomposes a unitary into:
I or T then a sequence of HT and/or SHT then a sequence of cliffords.

We use the VHeuristic (denominator exponent of DELTA) to guide the search.
VBHeuristic is only more effective in multi-qubit cases but MA_Form is only applicable for 1-qubit.
"""

RING = Ring.CYCLOTOMIC1
T = RING.get_gate("T-")
T_inv = RING.get_gate("T'-")
H = RING.get_gate("H-")
H_inv = RING.get_gate("H'-")
S = RING.get_gate("S-")
S_inv = RING.get_gate("S'-")
B = RING.get_gate("B-")
B_inv = RING.get_gate("B'-")

""""
Matusomoto-amano form requires three generators: T, HT, SHT
- we define them using gates from the CYCLOTOMIC ring.
"""

def get_matsumoto_gens(): #RETURNS list of generators, list of BW basis matrices

    gens = [
        Generator(T, T_inv,"T-", "T'-"),
        Generator(H*T, T_inv*H_inv, "HT-", "T'H'-"),
        Generator(S*H*T, T_inv*H_inv*S_inv, "SHT-", "T'H'S'-")]
    
    for g in gens:
        g.bmat = B_inv*g.mat_inv*B
    basis = [B_inv, B]

    return gens, basis


class MASearch(SearchClass):

    def __init__(self):
                
        cost_set, basis = get_matsumoto_gens()
        self.basis = basis
        self.gens = cost_set


    def go(self, start):

        #returns string sequence of generators and Clifford matrix 

        e_node = Vertex(self.basis[0]*start*self.basis[1])  #start node if start with I
        e_node.h = e_node.mat.denom_exp

        if e_node.h == 0:
            return e_node.sequence, start
        
        t_start_mat = self.gens[0].bmat * e_node.mat #start node if start with T = T_INV * START
        t_node = Vertex(t_start_mat)
        t_node.h = t_start_mat.denom_exp
        t_node.sequence = e_node.sequence + self.gens[0].symb 

        #Set the start either I or T
        if t_node.h<e_node.h:
            current_node = t_node
        else:
            current_node = e_node

        while current_node.h >0: 

            #expand current node by pre-multiplying by HT_inv, SHT_inv
            #and picking the one which reduces the denominator exponent (more)

            new_ht_mat = self.gens[1].bmat * current_node.mat
            new_sht_mat = self.gens[2].bmat * current_node.mat

            if new_ht_mat.denom_exp < new_sht_mat.denom_exp:
                new_node = Vertex(new_ht_mat)
                new_node.h = new_ht_mat.denom_exp
                new_node.sequence = current_node.sequence + self.gens[1].symb
                current_node = new_node
                print("HT gate added")

            else:
                new_node = Vertex(new_sht_mat)
                new_node.h = new_sht_mat.denom_exp
                new_node.sequence = current_node.sequence + self.gens[2].symb
                current_node = new_node
                print("SHT gate added")

        return current_node.sequence, self.basis[1]*current_node.mat*self.basis[0]

    def check_sequence(self,V,sequence,mat): 
            #overridden from SearchClass since MA_Form remainder clifford is at end of sequence, not start
    
            #gens_dict must be a dictionary of generators with g.symb as key
            #search through gens_dict for key matching start of sequence and post-multiply to recover U
            display(V)
            display(sequence)
            while sequence != "I-":
                if sequence.endswith("SHT-"):
                    sequence = sequence[:-4]
                    gen = self.gens[2]
                elif sequence.endswith("HT-"):
                    sequence = sequence[:-3]
                    gen = self.gens[1]
                elif sequence.endswith("T-"):
                    sequence = sequence[:-2]
                    gen = self.gens[0]
                V=gen.mat *V
                display(sequence)
            display(V)
            return V == mat

    