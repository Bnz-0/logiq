import math
from logiq import *


rd=1/math.sqrt(2)

class sdc:
    #SuperDense Coding

    basis = Basis([
            [rd,0,0,rd],
            [0,rd*1j, rd*1j, 0],
            [0, rd*1j, -rd*1j, 0],
            [rd,0,0,-rd]
        ], '0231')


    def __init__(self, qbits):
        self.qs = qbits


    @staticmethod
    def init(n):
        n += 2 if n%2==0 else 1 #for padding!
        qs_a=[] ; qs_b=[]
        for _ in range(0,n,2):
            qa = qbit(1,0) ; qb = qbit(1,0)
            Op.H | qa
            q = qa @ qb
            Op.cnot | q
            qs_a.append(qa)
            qs_b.append(qb)
        return sdc(qs_a), sdc(qs_b)

    
    def prepare(self, bits): #prepare the qubits to send the bits
        # add padding
        bits.append(1)
        if (len(bits)-1)%2 == 0: bits.append(0)

        for i in range(0,len(bits), 2):
            if bits[i] == 0:
                if bits[i+1] == 0: #00
                    pass
                else: #01
                    Op.Z | self.qs[i//2]
            else:
                if bits[i+1] == 0: #10
                    Op.X | self.qs[i//2]
                else: #11
                    Op.Y*1j | self.qs[i//2]


    def receive(self, qs_a):
        bits=[]
        for i in range(len(qs_a)):
            m = (qs_a[i] @ self.qs[i]).measure(sdc.basis)
            if m == sdc.basis.ew[0]: bits+=[0,0]
            elif m == sdc.basis.ew[1]: bits+=[1,0]
            elif m == sdc.basis.ew[2]: bits+=[1,1]
            else: bits+=[0,1]
        
        # remove padding
        while bits.pop() == 0: pass

        return bits


    def __getitem__(self,i):
        return self.qs[i]

    def __len__(self):
        return len(self.qs)

    def __repr__(self):
        return str(self.qs)
