from enum import Enum

from graphproblem.clifford import Clifford
from graphproblem.generator import Generator
from graphproblem.gatebuilder import get_gates_for_ring

"""Ring enum defines the types of matrices we can synthesise
 -- each enum calls the relevant gate generator function from gatebuilder.py so that
  we have the necessary gates for synthesising matrices in that ring. 
  RING.build_problem_gates also defines the sets needed to generate the generators."""

"""Example:
RING = Ring.CYCLOTOMIC2
free_set, cost_set, basis = RING.build_problem_gates
T = RING.get_gate("T-")
"""

class Ring(Enum): 
    CYCLOTOMIC1 = (1, "Z8","1 qubit clifford + t") #1 qubit, entries in ring Z[w, 1/2] 
    #AStar
    
    CYCLOTOMIC2 = (2,"Z8","2 qubit clifford + t") #2 qubits, entries in ring Z[w, 1/2] 
    #AStar
    
    CYCLOTOMIC2_with_CT = (2,"Z8","2 qubit clifford + t with CT") #2 qubits, entries in ring Z[w, 1/2]
    #AStar, greedy

    GAUSSIAN2 = (2, "ZI","2 qubits rational clifford")#2 qubits, Z[i]
    #AStar, greedy
    
    ZSQRT2_2 = (2, "ZSQRT2","2 qubits real clifford") #2 qubits, Z[sqrt2] 
    #AStar, greedy,

    DYADIC3 = (3,"DYADIC",6,"2 qubits real, rational clifford") #3 qubits, Z[1/2]

    def __init__(self,num_qubits,ring,description):
        self.num_qubits = num_qubits
        self.ring = ring
        self.description = description #necessary of cyclotomic2 and cyclotomic2_with_ct are defined identically.

    @property
    def ring_gates(self): #calls function from buildgates.py and returns the dictionary of gates
        return get_gates_for_ring(self.ring, self.num_qubits)

    def get_gate(self, name): #returns the specific gate implemented in that ring
        return self.ring_gates[name] #naming format: "H-" calls Hadamard gate

    def build_problem_gates(self):
        #returns [free_set, cost_set, transformation matrices for Barnes-Wall basis]
        
        gates = self.ring_gates
        EYE = gates["I-"]
        B = gates["B-"]
        B_inv = gates["B'-"]
        
        if self is Ring.DYADIC3:
            X = gates["X-"]
            Z = gates["Z-"]

            CX = gates["CX-"]
            XC = gates["XC-"]
            CZ = gates["CZ-"]
            ZC = gates["CZ-"]

            IHH = gates["IHH-"]
            HHI = gates["HHI-"]
            HIH = gates["HIH-"]

            #IXI = EYE@X@EYE
            #XII = X@EYE@EYE
            IIX = EYE@EYE@X
            ZII = Z@EYE@EYE
            IZI = EYE@Z@EYE
            IIZ= EYE@EYE@Z

            ICX = EYE@CX
            IXC = EYE@XC
            #CXI = CX@EYE
            XCI  = XC@EYE
            #XIC = gates["XIC-"]
            #CIX = gates["CIX-"]

            ICZ = gates["ICZ-"]
            #IZC = EYE@ZC
            CZI = gates["CZI-"]
            #ZCI = ZC@EYE
            CIZ = gates["CIZ-"]
            #ZIC = gates["ZIC-"]

            #CCX= gates["CCX-"]
            CCZ = gates["CCZ-"]

            B2 = gates["B2-"]
            B2_inv = gates["B2'-"]

            return ( 
                [
                    Clifford(HHI, HHI, "HHI-", "HHI'-"),
                    Clifford(HIH, HIH, "HIH-", "HIH'-"),
                    #Clifford(IHH, IHH, "IHH-", "IHH'-"), #dont need all three, since HHI * HIH = IHH

                    #Clifford(XII, XII, "XII-", "XII'-"),
                    #Clifford(IXI, IXI, "IXI-", "IXI'-"),
                    #Clifford(IIX, IIX, "IIX-", "IIX'-"),

                    Clifford(ZII, ZII, "ZII-", "ZII'-"),
                    Clifford(IZI, IZI, "IZI-", "IZI'-"),
                    Clifford(IIZ, IIZ, "IIZ-", "IIZ'-"),

                    Clifford(CZI, CZI, "CZI-", "CZI'-"),
                    Clifford(ICZ, ICZ, "ICZ-", "ICZ'-"),
                    Clifford(CIZ, CIZ, "CIZ-", "CIZ'-"),

                    #Clifford(XCI, XCI, "XCI-", "XCI'-"),
                    #Clifford(IXC, IXC, "IXC-", "IXC'-"),
                    #Clifford(XIC, XIC, "XIC-", "XIC'-"),
                    #Clifford(CXI, CXI, "CXI-", "CXI'-"),
                    #Clifford(ICX, ICX, "ICX-", "ICX'-"),
                    #Clifford(CIX, CIX, "CIX-", "ZIX'-"),
                ],
                [
                    #Generator(CCX, CCX, "CCX-", "CCX'-"),
                    #Generator(CXC, CXC, "CXC-", "CXC'-"),
                    #Generator(XCC, XCC, "XCC-", "XCC'-"),
                    Generator(CCZ, CCZ, "CCZ-", "CCZ'-"),
                    #Generator(CZC, CZC, "CZC-", "CZC'-"),
                    #Generator(ZCC, ZCC, "ZCC-", "ZCC'-"),
                ],
                [B_inv, B, B2_inv, B2]) #first basis CCZ is unimodular, second basis is CCX unimodular
            
        #gates with an implementation in every following ring 
        H = gates["H-"]
        H_inv = gates["H'-"]

        if self.ring == "Z8":
            S = gates["S-"]
            S_inv = gates["S'-"]
            T = gates["T-"]
            T_inv = gates["T'-"]

            if self is Ring.CYCLOTOMIC1:
                return (
                    [Clifford(H, H_inv, "H-", "H'-"), 
                     Clifford(S, S_inv, "S-", "S'-")
                     ],
                    [Generator(T, T_inv, "T-", "T'-")
                     ],
                    [B_inv, B],
                )

            CX = gates["CX-"]
            CX_inv = gates["CX'-"]

            IH = EYE @ H
            HI = H @ EYE
            IT = EYE @ T
            TI = T @ EYE
            SI = S @ EYE
            IS = EYE @ S
    
            SI_inv = S_inv @ EYE
            IS_inv = EYE @ S_inv
            TI_inv = T_inv @ EYE
            IT_inv = EYE @ T_inv
            HI_inv = H_inv@ EYE
            IH_inv = EYE@H_inv


            if self is Ring.CYCLOTOMIC2:
                return ( 
                [
                    Clifford(IH, IH_inv, "IH-", "I'H'-"),
                    Clifford(HI, HI_inv, "HI-", "H'I'-"),
                    Clifford(SI, SI_inv, "SI-", "S'I'-"),
                    Clifford(IS, IS_inv, "IS-", "I'S'-"),
                    Clifford(CX, CX_inv, "CX-", "CX'-"),
                ],
                [
                    Generator(IT, IT_inv, "IT-", "IT'-"),
                    Generator(TI, TI_inv, "TI-", "TI'-")
                ],
                [B_inv @ B_inv, B @ B])

            elif self is Ring.CYCLOTOMIC2_with_CT:

                CT = gates["CT-"]
                CT_inv =gates["CT'-"]
                return ( 
                    [
                        Clifford(IH, IH_inv, "IH-", "I'H'-"),
                        Clifford(HI, HI_inv, "HI-", "H'I'-"),
                        Clifford(SI, SI_inv, "SI-", "S'I'-"),
                        Clifford(IS, IS_inv, "IS-", "I'S'-"),
                        Clifford(CX, CX_inv, "CX-", "CX'-")
                    ],
                    [
                        Generator(IT, IT_inv, "IT-", "IT'-"),
                        Generator(TI, TI_inv, "TI-", "TI'-"),
                        Generator(CT,CT_inv, "CT-", "CT'-")
                    ],
                    [B_inv @ B_inv, B @ B])



        elif self is Ring.GAUSSIAN2:
            S = gates["S-"]
            S_inv = gates["S'-"]
            CS = gates["CS-"]
            CS_inv = gates["CS'-"]
            CX = gates["CX-"]
            CX_inv = gates["CX'-"]
    
            IH = EYE @ H
            HI = H @ EYE
            SI = S @ EYE
            IS = EYE @ S
    
            SI_inv = S_inv @ EYE
            IS_inv = EYE @ S_inv
            HI_inv = H_inv@ EYE
            IH_inv = EYE@H_inv

            return (
            [
                Clifford(IH, IH_inv, "IH-", "I'H'-"),
                Clifford(HI, HI_inv, "HI-", "H'I'-"),
                Clifford(SI, SI_inv, "SI-", "S'I'-"),
                Clifford(IS, IS_inv, "IS-", "I'S'-"),
                Clifford(CX, CX_inv, "CX-", "CX'-")
            ],
            [
                Generator(CS, CS_inv, "CS-", "CS'-"),
            ],
            [B_inv @ B_inv, B @ B])

        elif self is Ring.ZSQRT2_2:
            X = gates["X-"]
            Z = gates["Z-"]
            X_inv = gates["X'-"]
            Z_inv = gates["Z'-"]
            CX = gates["CX-"]
            CX_inv = gates["CX'-"]
            CH = gates["CH-"]
            CH_inv = gates["CH'-"]
            HC = gates["HC-"]
            HC_inv = gates["HC'-"]

            IH = EYE @ H
            HI = H @ EYE
            XI = X @ EYE
            IX = EYE @ X
            ZI = Z @ EYE
            IZ = EYE @ Z 
    
            XI_inv = X_inv @ EYE
            IX_inv = EYE @ X_inv
            HI_inv = H_inv@ EYE
            IH_inv = EYE@H_inv
            ZI_inv = Z_inv@ EYE
            IZ_inv = EYE@Z_inv

            return (
                [
                    Clifford(IH, IH_inv, "IH-", "I'H'-"),
                    Clifford(HI, HI_inv, "HI-", "H'I'-"),
                    Clifford(XI, XI_inv, "SI-", "S'I'-"),
                    Clifford(IX, IX_inv, "IS-", "I'S'-"),
                    Clifford(ZI, ZI_inv, "ZI-", "ZI'-"),
                    Clifford(IZ, IZ_inv, "IZ-", "IZ'-"),
                    Clifford(CX, CX_inv, "CX-", "CX'-"),
                ],
                [
                    Generator(CH, CH_inv, "CH-", "CH'-"),
                    Generator(HC, HC_inv, "HC-", "HC'-")
                ],
                [B_inv @ B_inv, B @ B])
        else:
            raise ValueError(f"No gate set definition defined for {self.name}")
        
