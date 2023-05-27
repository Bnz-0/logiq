from .abs_Qstate import _Qstate, need_norm, unreal
from .Basis import Basis, CanonBasis
from .Operator import MeasureOp, Op
from .Qbits import Qbits
from .Qerrors import IllegalOperationError, InitializationError
from .Qmath import roundedVector, vector
from .qtils import equal, formatProbs, val2str


#### Qbit.py
#
# This file contains the Qbit and Qbit_ent classes.
# This 2 classes serve to create a pure quantum state even in an entangled state,
# i.e. a single qubit (or better: every kind of qudit)
#
####


# ↓↓↓↓↓↓↓↓↓↓↓↓ Qbit class ↓↓↓↓↓↓↓↓↓↓↓↓ #

class Qbit(_Qstate):
    
    # Single qubit class

    def __init__(self, state, basis = None, normalize = None, transform = False):
        try:
            super().__init__(vector(state), basis, [self], len(state))
            self._ent = None

            if not equal(self._state.norm(), 1):
                if (normalize is None and need_norm()) or normalize:
                    self._state.normalize()
                else:
                    raise IllegalOperationError('A quantum state must have norm 1, this one have norm '+str(self._state.norm()))
            
            if self._basis is None: self.setBasis(None)
            if transform: self._state = self._basis.transform(self._state)

        except Exception as e:
            raise InitializationError("Error to initialize Qbit", e)



    def setBasis(self, basis):
        if basis is None:
            self._basis = CanonBasis(len(self))
        else:
            self._basis = basis
            if len(self) != len(self._basis):
                raise ValueError('Lengths of state and basis must be equal')


    def _entangle(self, ent, pos):
        Qbit_ent.__init__(self, ent, pos, self._basis)
        self.__class__ = Qbit_ent


    def measure(self, basis = None):
        i, _, basis = super().measure(basis)
        self._state = vector(basis[i])

        return basis.ew[i]


    def __matmul__(self, q):
        if isinstance(q, int):
            return Qbits([Qbit(vector(self._state), self._basis) for _ in range(q)])
        return Qbits([self, q])


    @unreal
    def printAs(self, basis):
        s=''
        state = roundedVector(basis.transform(self._getState()))
        for i in range(len(self)):
            if(not equal(state[i], 0)):
                s += val2str(state[i]) + '|' + basis.symbols[i] + '> '

        return s[:-1]


    @unreal
    def printProbs(self, basis = None):
        b = basis if basis is not None else self._basis
        print(formatProbs(b.transform(self._getState()), b.symbols))

# ↑↑↑↑↑↑↑↑↑↑↑↑ Qbit class ↑↑↑↑↑↑↑↑↑↑↑↑ #



# ↓↓↓↓↓↓↓↓↓↓↓↓ Qbit_ent class ↓↓↓↓↓↓↓↓↓↓↓↓ #

class Qbit_ent(Qbit, _Qstate):

    # Single qubit class for entangled states

    def __init__(self, ent, pos, basis):
        _Qstate.__init__(self, None, basis, [self], self._length)
        self._ent = ent
        self._pos = pos


    def _getState(self):
        return self._ent._calc_p_state(self._pos)


    def _entangle(self, ent, pos):
        self._ent = ent
        self._pos = pos


    def apply(self, op, pos = None):
        self._ent.apply(op, self._pos)


    def measure(self, basis = None):
        i, state, basis = _Qstate.measure(self, basis)
        self._ent.apply(MeasureOp(state, basis, i), self._pos)

        return basis.ew[i]



    def __matmul__(self, q):
        if isinstance(q, int):
            raise IllegalOperationError('This state is in entanglement')
        return Qbits([self, q])
    
    
    def printAs(self, basis):
        s='' ; rappr = '|>'
        state = roundedVector(basis.transform(self._getState()))
        for i in range(len(self)):
            if(not equal(state[i], 0)):
                s += val2str(state[i]) + rappr[0] + basis.symbols[i] + rappr[1] + ' '
        
        s = s.replace('+', '±').replace('-', '±')
        return s[:-1]

# ↑↑↑↑↑↑↑↑↑↑↑↑ Qbit_ent class ↑↑↑↑↑↑↑↑↑↑↑↑ #
