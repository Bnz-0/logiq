import math
from logiq import *


def deutch(f):
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



def grover(db, f):
    vect = lambda x: ket(*(1 if i==x else 0 for i in range(len(db))))
    Uw = Op.build(
        {i: -vect(i) if f(i) else vect(i) for i in range(len(db))}
    )
    Us = Op(2*db.state * ~db.state - matrix.Id(len(db)))

    G = Uw*Us

    for i in range(int(round(math.sqrt(len(db))))):
        G | db

    db.measure(CanonBasis(len(db)))
    return db
