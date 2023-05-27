from random import random

from .Basis import Basis
from .Operator import Op
from .Qerrors import IllegalOperationError, NotAllowError
from .Qmath import mod_square, vector


#### abs_Qstate.py
#
# This file contains the abstract class for Qstate
# and the functions to manage the auto-normalization and the cheating permissions
#
####


# ↓↓↓↓↓↓↓↓↓↓↓↓ Auto-normalization and cheating permissions functions ↓↓↓↓↓↓↓↓↓↓↓↓ #

_auto_norm = False
_cheat = True


def auto_norm(value = None):
    "If `True`, when a qubit is created it will be normalize if necessary"
    global _auto_norm
    _auto_norm = value if value in (False, True) else not _auto_norm


def cheat(value):
    "If `False`, disables the ability to see what is inside a qubit"
    if value in (False, True):
        global _cheat
        _cheat = value


def can_cheat():
    return _cheat


def need_norm():
    return _auto_norm


def unreal(func):
    # according to _cheat value, allow or not to launch the function that it decorates
    def wrapper(*args, **kwargs):
        if not _cheat:
            raise NotAllowError("Cheating is not allowed!")
        return func(*args, **kwargs)
    return wrapper

# ↑↑↑↑↑↑↑↑↑↑↑↑ Auto-normalization and cheating permissions functions ↑↑↑↑↑↑↑↑↑↑↑↑ #



# ↓↓↓↓↓↓↓↓↓↓↓↓ _Qstate class ↓↓↓↓↓↓↓↓↓↓↓↓ #


class _Qstate:
    # Abstract class for Qstate

    def __init__(self, state, basis, qbits, length):
        self._state = state
        self._basis = basis
        self._qbits = qbits
        self._length = length
        
    
    @property
    @unreal
    def state(self):
        return self._getState()
    

    @property
    def basis(self):
        return self._basis


    @property
    def spaces(self):
        "Return the number of qubits that compose this quantum state (i.e. the number of vectorial spaces)"
        return len(self._qbits)


    def __getitem__(self, i):
        return self._qbits[i]


    def _getState(self):
        return self._state


    def setBasis(self, basis):
        "Set the default basis for this quantum state"
        pass
    

    def apply(self, op, pos=None):
        "Apply an operator to this quantum state"
        self._state = op * self._state
        self._state.round_error((0,1))


    def apply2all(self, op):
        "Apply the operator to all the qubit that compose this quantum state"
        for q in self._qbits:
            q.apply(op)


    def measure(self, basis = None):
        "Measure this state in Basis `basis`, the result will be the eigenvalue associated"
        if basis is None:
            if self._basis is None: raise IllegalOperationError('Measurement basis required')
            basis = self._basis
        state = basis.transform(self._getState())
        r = random()
        i = -1 ; s = 1
        while s > r:
            i += 1
            s -= mod_square(state[i])
        return i, state, basis


    @unreal
    def prob(self, i, basis = None):
        "Returns the probability to measure the `i`-th state if the basis `basis` is used"
        if basis is None: basis = self._basis
        if basis is None: raise IllegalOperationError('Need a basis to calculate the probabilities')
        M = (basis[i] * ~basis[i])
        return ~self._getState() * (~M * M) * self._getState()


    @unreal
    def printProbs(self, basis = None):
        "Returns a string that describe the probabilities of measurement using basis `basis` of this quantum state"
        pass


    @unreal
    def printAs(self, basis):
        "Returns the representation of this quantum state from the 'point of view' of `basis`"
        pass


    def __len__(self):
        return self._length


    def __hash__(self):
        return id(self)
    

    def __str__(self):
        return self.printAs(self._basis)


    def __repr__(self):
        return self.__str__() if can_cheat() else '<Qstate @ '+hex(id(self))+', ('+self.__class__.__qualname__+')>'

# ↑↑↑↑↑↑↑↑↑↑↑↑ _Qstate class ↑↑↑↑↑↑↑↑↑↑↑↑ #
