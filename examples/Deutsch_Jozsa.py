#!/bin/env python3
from logiq import Op, Basis, qbit

# https://en.wikipedia.org/wiki/Deutsch%E2%80%93Jozsa_algorithm

def Deutsch_Jozsa(f):
    Uf = Op.build({			#          f(0)=f(1)=0   f(0)=f(1)=1  f(0)=0,f(1)=1  f(0)=1,f(1)=0
            0: 0^f(0),		# |00>  -->    |00>          |01>          |00>          |01>
            1: 1^f(0),		# |01>  -->    |01>          |10>          |01>          |00>
            2: 2+0^f(1),	# |10>  -->    |10>          |11>          |11>          |10>
            3: 2+1^f(1)		# |11>  -->    |11>          |10>          |10>          |11>
        })
    B = Basis([[1,0],[0,-1]])
    q = qbit(1,0, basis=B) @ 2 #2 qubits having the same state: |00>

    Op.X | q[1]
    Op.H ^ q
    Uf | q
    Op.H | q[0]

    return q[0].measure(B)

def main():
    r = Deutsch_Jozsa(lambda x: x)
    print(f"Deutsch_Jozsa(f(x) = x) --> {r}")
    assert(r == -1)
    r = Deutsch_Jozsa(lambda x: 0)
    print(f"Deutsch_Jozsa(f(x) = 0) --> {r}")
    assert(r == 1)

if __name__ == "__main__":
    main()
