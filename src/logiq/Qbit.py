from .abs_Qstate import _Qstate, need_norm, unreal
from .Basis import CanonBasis
from .Operator import MeasureOp
from .Qbits import Qbits
from .Qerrors import IllegalOperationError, InitializationError
from .Qmath import rounded_vector, Vector
from .qtils import equal, format_probs, val2str


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
            super().__init__(Vector(state), basis, [self], len(state))
            self._ent = None

            if not equal(self._state.norm(), 1):
                if (normalize is None and need_norm()) or normalize:
                    self._state.normalize()
                else:
                    raise IllegalOperationError('A quantum state must have norm 1, this one have norm '+str(self._state.norm()))

            if self._basis is None: self.set_basis(None)
            if transform: self._state = self._basis.transform(self._state)

        except Exception as e:
            raise InitializationError("Error to initialize Qbit", e) from e

    def set_basis(self, basis=None):
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
        self._state = Vector(basis[i])
        return basis.ew[i]

    def __matmul__(self, q):
        if isinstance(q, int):
            return Qbits([Qbit(Vector(self._state), self._basis) for _ in range(q)])
        return Qbits([self, q])

    @unreal
    def print_as(self, basis):
        s = ''
        state = rounded_vector(basis.transform(self._get_state()))
        for i in range(len(self)):
            if not equal(state[i], 0):
                s += val2str(state[i]) + '|' + basis.symbols[i] + '> '

        return s[:-1]

    @unreal
    def print_probs(self, basis = None):
        b = basis if basis is not None else self._basis
        print(format_probs(b.transform(self._get_state()), b.symbols))

# ↑↑↑↑↑↑↑↑↑↑↑↑ Qbit class ↑↑↑↑↑↑↑↑↑↑↑↑ #


# ↓↓↓↓↓↓↓↓↓↓↓↓ Qbit_ent class ↓↓↓↓↓↓↓↓↓↓↓↓ #

class Qbit_ent(Qbit, _Qstate):

    # Single qubit class for entangled states

    def __init__(self, ent, pos, basis):
        _Qstate.__init__(self, None, basis, [self], self._length)
        self._ent = ent
        self._pos = pos

    def _get_state(self):
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

    def print_as(self, basis):
        s = ''
        rappr = '|>'
        state = rounded_vector(basis.transform(self._get_state()))
        for i in range(len(self)):
            if not equal(state[i], 0):
                s += val2str(state[i]) + rappr[0] + basis.symbols[i] + rappr[1] + ' '

        return s.replace('+', '±').replace('-', '±')[:-1]

# ↑↑↑↑↑↑↑↑↑↑↑↑ Qbit_ent class ↑↑↑↑↑↑↑↑↑↑↑↑ #
