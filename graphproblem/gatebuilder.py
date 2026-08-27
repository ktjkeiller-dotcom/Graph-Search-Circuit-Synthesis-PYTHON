from matrixnormalforms.eds import z8, zi, zsqrt2, z
from matrixnormalforms import q
from matrixnormalforms.qmatrix import QMatrix as fm

"""This holds the constructions of necessary gates to synthesise matrices over a given ring.

When setting up a RING enum from rings.py, we call RING.build_problem_gates, which returns the list of 
free_gates, cost_gates, and basis_change gates.

Functions "build_z..._gates():"
- must define the implementations of basic integers in the ring, eg. 1, 0, w, i, sqrt2
- this then allows for general constructions of gates shared between rings

Note:
Some gates have different constructions in different rings. 
- The complex rings have a complex Hadamard gate (phase shifted from the standard) 
but the real, irrational rings have the standard Hadamard gate
- The BW basis change matrices are different between rings depending on if they are
complex or irrational
 """

_RING_BUILDERS = {}

def make_gates(vals, ring, qubits): #returns a dictionary of basic gates that can be implemented in the ring. 
    gates = {}

    #gates that every set implements in the same way
    gates["I-"] = fm(2,2,[vals["1"],vals["0"],vals["0"],vals["1"]])
    gates["X-"] = fm(2,2,[vals["0"],vals["1"],vals["1"],vals["0"]])
    gates["X'-"] = gates["X-"].copy()
    gates["Z-"] = fm(2,2,[vals["1"],vals["0"],vals["0"],-vals["1"]])
    gates["Z'-"] = gates["Z-"].copy()

    if qubits>1:
        gates["CX-"] = fm(4,4,[vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["1"],vals["0"]])
        gates["CX'-"] = gates["CX-"].copy()
        gates["CZ-"] = fm(4,4,[vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],-vals["1"]])
        gates["CZ'-"] = gates["CZ-"].copy()
        gates["XC-"] = fm(4,4,[vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"]])
        gates["XC'-"] = gates["XC-"].copy()

        swap = fm(4,4,
            [vals["1"], vals["0"], vals["0"], vals["0"],
            vals["0"], vals["0"], vals["1"], vals["0"],
            vals["0"], vals["1"], vals["0"], vals["0"],
            vals["0"], vals["0"], vals["0"], vals["1"]])
        
    if "i" in vals: #define complex gates
        gates["S-"] = fm(2,2,[vals["1"],vals["0"],vals["0"],vals["i"]])
        gates["S'-"] = fm(2,2,[vals["1"],vals["0"],vals["0"],-vals["i"]])
        gates["H-"] = fm(2,2,[vals["1/1+i"],vals["1/1+i"],vals["1/1+i"],-vals["1/1+i"]])
        gates["H'-"] = fm(2,2,[vals["1/1-i"],vals["1/1-i"],vals["1/1-i"],-vals["1/1-i"]])
        gates["B-"] = fm(2,2,[vals["1/1+i"],vals["0"],vals["1/1+i"],vals["1"]])
        gates["B'-"] = fm(2,2,[vals["1"]+vals["i"],vals["0"],-vals["1"],vals["1"]])

        if qubits >1:
            gates["CS-"] = fm(4,4,[vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["i"]])
            gates["CS'-"] = fm(4,4,[vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],-vals["i"]])
            gates["CH-"] = fm(4,4,[vals["1"],vals["0"], vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1/1+i"],vals["1/1+i"],vals["0"],vals["0"],vals["1/1+i"],-vals["1/1+i"]])
            gates["CH'-"] = fm(4,4,[vals["1"],vals["0"], vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1/1-i"],vals["1/1-i"],vals["0"],vals["0"],vals["1/1-i"],-vals["1/1-i"]])       
                       

    if ring == "ZSQRT2": #define gates specific or needed in this ring
        gates["H-"] = fm(2,2,[vals["1/sqrt2"],vals["1/sqrt2"],vals["1/sqrt2"],-vals["1/sqrt2"]])
        gates["H'-"] = fm(2,2,[vals["1/sqrt2"],vals["1/sqrt2"],vals["1/sqrt2"],-vals["1/sqrt2"]])
        gates["B-"] = fm(2,2,[vals["1/sqrt2"],vals["0"],vals["1/sqrt2"],vals["1"]])
        gates["B'-"] = fm(2,2,[vals["sqrt2"],vals["0"],-vals["1"],vals["1"]])

        if qubits >1:
            gates["CH-"] = fm(4,4,[vals["1"],vals["0"], vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1/sqrt2"],vals["1/sqrt2"],vals["0"],vals["0"],vals["1/sqrt2"],-vals["1/sqrt2"]])
            gates["CH'-"] = gates["CH-"].copy()
            gates["HC-"] = swap * gates["CH-"] * swap
            gates["HC'-"] = swap * gates["CH'-"] * swap
                    

    if ring == "Z8": #define gates specific to cyclotomic set
        gates["T-"] = fm(2,2,[vals["1"],vals["0"],vals["0"],vals["w"]])
        gates["T'-"] = fm(2,2,[vals["1"],vals["0"],vals["0"],-vals["w"]*vals["i"]])

        if qubits>1:
            gates["CT-"] = fm(4,4,[vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["w"]])
            gates["CT'-"] = fm(4,4,[vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],-vals["w"]*vals["i"]])
            

    if ring == "DYADIC": #define gates specific to dyadic ring 
        gates["B-"] = fm(8,8,[
            vals["1/2"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],
            vals["0"],vals["1/2"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],
            vals["0"],vals["0"],vals["1/2"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],
            vals["1/2"],vals["1/2"],vals["1/2"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],
            vals["0"],vals["0"],vals["0"],vals["0"],vals["1/2"],vals["0"],vals["0"],vals["0"],
            vals["1/2"],vals["1/2"],vals["0"],vals["0"],vals["1/2"],vals["1"],vals["0"],vals["0"],
            vals["1/2"],vals["0"],vals["1/2"],vals["0"],vals["1/2"],vals["0"],vals["1"],vals["0"],
            vals["0"],vals["1/2"],vals["1/2"],vals["0"],vals["1/2"],vals["0"],vals["0"],vals["1"]
        ])

        gates["B'-"] = fm(8,8,[
            vals["2"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],
            vals["0"],vals["2"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],
            vals["0"],vals["0"],vals["2"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],
            -vals["1"],-vals["1"],-vals["1"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],
            vals["0"],vals["0"],vals["0"],vals["0"],vals["2"],vals["0"],vals["0"],vals["0"],
            -vals["1"],-vals["1"],vals["0"],vals["0"],-vals["1"],vals["1"],vals["0"],vals["0"],
            -vals["1"],vals["0"],-vals["1"],vals["0"],-vals["1"],vals["0"],vals["1"],vals["0"],
            vals["0"],-vals["1"],-vals["1"],vals["0"],-vals["1"],vals["0"],vals["0"],vals["1"]
        ])

        gates["B2-"] = fm(8,8,[
            vals["1/2"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],
            vals["1/2"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],
            vals["1/2"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],
            vals["1/2"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],
            vals["1/2"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],
            vals["1/2"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],
            vals["1/2"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],
            vals["1/2"],vals["1"],vals["1"],vals["1"],vals["1"],vals["1"],vals["1"],vals["2"]
        ])

        gates["B2'-"] = fm(8,8,[
            vals["2"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],
            -vals["1"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],
            -vals["1"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],
            -vals["1"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],
            -vals["1"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],vals["0"],
            -vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],vals["0"],
            -vals["1"],vals["0"],vals["0"],vals["0"],vals["0"],vals["0"],vals["1"],vals["0"],
            vals["5/2"],-vals["1/2"],-vals["1/2"],-vals["1/2"],-vals["1/2"],-vals["1/2"],-vals["1/2"],vals["1/2"]
            ])

        gates["HH-"] = fm(4,4,[
            vals["1/2"],  vals["1/2"],  vals["1/2"], vals["1/2"],
            vals["1/2"], -vals["1/2"],  vals["1/2"],-vals["1/2"],
            vals["1/2"],  vals["1/2"], -vals["1/2"],-vals["1/2"],
            vals["1/2"], -vals["1/2"], -vals["1/2"], vals["1/2"]
        ])
        gates["HH'-"] = gates["HH-"].copy()

        gates["IHH-"] = gates["I-"] @ gates["HH-"]
        swap_1_2 = swap @ gates["I-"]
        gates["HIH-"] = swap_1_2 * gates["IHH-"] * swap_1_2
        gates["HHI-"] = gates["HH-"] @ gates["I-"]

        gates["CZI-"] = gates["CZ-"] @ gates["I-"]
        gates["ICZ-"] = gates["I-"] @ gates["CZ-"]
        gates["CIZ-"] = swap_1_2 * gates["ICZ-"] * swap_1_2

        gates["CCX-"] = fm(8,8,[
            vals["1"], vals["0"], vals["0"], vals["0"],vals["0"], vals["0"], vals["0"], vals["0"],
            vals["0"], vals["1"], vals["0"], vals["0"],vals["0"], vals["0"], vals["0"], vals["0"],
            vals["0"], vals["0"], vals["1"], vals["0"],vals["0"], vals["0"], vals["0"], vals["0"],
            vals["0"], vals["0"], vals["0"], vals["1"],vals["0"], vals["0"], vals["0"], vals["0"],
            vals["0"], vals["0"], vals["0"], vals["0"],vals["1"], vals["0"], vals["0"], vals["0"],
            vals["0"], vals["0"], vals["0"], vals["0"],vals["0"], vals["1"], vals["0"], vals["0"],
            vals["0"], vals["0"], vals["0"], vals["0"],vals["0"], vals["0"], vals["0"], vals["1"],
            vals["0"], vals["0"], vals["0"], vals["0"],vals["0"], vals["0"], vals["1"], vals["0"]
            ])

        gates["CCZ-"] = fm(8,8,[
                vals["1"], vals["0"], vals["0"], vals["0"],vals["0"], vals["0"], vals["0"], vals["0"],
                vals["0"], vals["1"], vals["0"], vals["0"],vals["0"], vals["0"], vals["0"], vals["0"],    
                vals["0"], vals["0"], vals["1"], vals["0"],vals["0"], vals["0"], vals["0"], vals["0"],   
                vals["0"], vals["0"], vals["0"], vals["1"],vals["0"], vals["0"], vals["0"], vals["0"],    
                vals["0"], vals["0"], vals["0"], vals["0"],vals["1"], vals["0"], vals["0"], vals["0"],   
                vals["0"], vals["0"], vals["0"], vals["0"],vals["0"], vals["1"], vals["0"], vals["0"],   
                vals["0"], vals["0"], vals["0"], vals["0"],vals["0"], vals["0"], vals["1"], vals["0"],   
                vals["0"], vals["0"], vals["0"], vals["0"],vals["0"], vals["0"], vals["0"], -vals["1"]
                ])

    return gates 

def build_z8_gates(ring, qubits):
    vals = {}
    vals["1"] = q.Q([1,0,0,0],0,z8.Z8)
    vals["w"] = q.Q([0,1,0,0],0,z8.Z8)
    vals["i"] = q.Q([0,0,1,0],0,z8.Z8)
    vals["0"] = q.Q([0,0,0,0],0,z8.Z8)
    vals["sqrt2_inv"] = q.Q([1,1,1,0],2,z8.Z8) #1/sqrt2
    vals["sqrt2"] = q.Q([0,1,0,-1],0,z8.Z8)
    vals["1/1+i"] = q.Q([1,1,0,-1],2,z8.Z8) #1/(1+i)
    vals["1/1-i"] = q.Q([0,1,1,1],2,z8.Z8)

    gates = make_gates(vals, ring, qubits)
    return gates

def build_zi_gates(ring, qubits):
    vals = {}
    vals["1"] = q.Q([1,0],0,zi.ZI)
    vals["i"] = q.Q([0,1],0,zi.ZI)
    vals["0"] = q.Q([0,0],0,zi.ZI)
    vals["1/1+i"] = q.Q([1,0],1,zi.ZI) #1/(1+i)
    vals["1/1-i"] = q.Q([0,-1],1,zi.ZI) #1/(1-i)

    gates = make_gates(vals, ring, qubits)
    return gates

def build_zsqrt2_gates(ring, qubits):
    vals = {}
    vals["1"] = q.Q([1,0],0,zsqrt2.ZSQRT2)
    vals["sqrt2"] = q.Q([0,1],0,zsqrt2.ZSQRT2)
    vals["0"] = q.Q([0,0],0,zsqrt2.ZSQRT2)
    vals["1/sqrt2"] = q.Q([1,0],1,zsqrt2.ZSQRT2) #1/(1+i)
    vals["2"] =q.Q([2,0],0,zsqrt2.ZSQRT2)
    vals["1/2"] = q.Q([1,0],2, zsqrt2.ZSQRT2)

    gates = make_gates(vals, ring, qubits)
    return gates

def build_dyadic_gates(ring, qubits):
    vals = {}
    vals["1"] = q.Q([1],0,z.Z)
    vals["0"] = q.Q([0],0,z.Z)
    vals["1/2"] = q.Q([1],1,z.Z)
    vals["2"] = vals["1"]+vals["1"]
    vals["5/2"] = vals["2"] + vals["1/2"]

    gates = make_gates(vals, ring, qubits)
    return gates

_RING_BUILDERS["Z8"] = build_z8_gates
_RING_BUILDERS["ZI"] = build_zi_gates
_RING_BUILDERS["ZSQRT2"] = build_zsqrt2_gates
_RING_BUILDERS["DYADIC"] = build_dyadic_gates

def get_gates_for_ring(ring, qubits):
    if ring not in _RING_BUILDERS:
        raise KeyError(f"Unknown ring: {ring}")
    return _RING_BUILDERS[ring](ring, qubits)


def get_ring_gate(ring_name, gate_name):
    return get_gates_for_ring(ring_name)[gate_name]