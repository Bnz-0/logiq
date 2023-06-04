from .abs_Qstate import _Qstate, unreal
from .Operator import MeasureOp, Op
from .Qerrors import IllegalOperationError, InitializationError
from .Qmath import ket, math, Matrix, rounded_vector
from .qtils import Vdigit, equal, format_probs, mod_square, prod, val2str

#### Qbits.py
#
# This file contains the Qbits and Qent classes.
# These 2 classes serve to create composed quantum states,
# i.e. multiple qubit (or better: every kind of qudit) states, even entangled.
#
# The Qent classes is a particular class that manage (hiddenly) the entangled states
#
####


# ↓↓↓↓↓↓↓↓↓↓↓↓ Common function ↓↓↓↓↓↓↓↓↓↓↓↓ #

def complete_Op(qbits, operator, pos):
    # it builds the operator to be applied, according to the given operator and the start position
    real_op = Op.neutral()
    if pos is None: pos=0
    if isinstance(operator, (list, tuple)): #set of operators
        if isinstance(pos, (list, tuple)):
            for i in range(len(qbits)):
                real_op @= operator[i] if i in pos else Op.Id(len(qbits[i]))
        else:
            for i in range(pos): real_op @= Op.Id(len(qbits[i]))
            for op in operator: real_op @= op

    else: # single operators
        i=0
        while i<len(qbits):
            if i==pos:
                real_op @= operator
                l=1
                while l<len(operator):
                    l*=len(qbits[i])
                    i+=1
            else:
                real_op @= Op.Id(len(qbits[i]))
                i+=1

    return real_op


def is_ent(q):
    return q._ent is not None


def get_ent(qents):
    # get the Qbits_ent common to all Qbit in qents
    ent = qents[0]._ent
    if ent is None: return None
    for q in qents:
        if ent is not q._ent:
            ent=None
            break

    return ent


def gen_ent(qents):
    # if doesn't exist a common Qbits_ent, this generates it
    qent = get_ent(qents)
    if qent is None:
        qent = Qent(qents)
    return qent

# ↑↑↑↑↑↑↑↑↑↑↑↑ Common function ↑↑↑↑↑↑↑↑↑↑↑↑ #


# ↓↓↓↓↓↓↓↓↓↓↓↓ Qbits class ↓↓↓↓↓↓↓↓↓↓↓↓ #

class Qbits(_Qstate):

    # Multiple qubit class

    def __init__(self, qbits, basis=None):
        try:
            qs = []
            for qbit in qbits:
                if isinstance(qbit, Qbits):
                    qs.extend(qbit._qbits)
                else:
                    qs.append(qbit)

            if len(qs) != len(set(qs)):
                raise IllegalOperationError('Duplicate qubit not allowed')

            super().__init__(None, basis, qs, prod((len(q) for q in qs)))

        except Exception as e:
            raise InitializationError('Error to initialize the Qbits', e) from e

    def set_basis(self, basis = None):
        self._basis = basis

    def _get_state(self):
        ent = get_ent(self._qbits)
        # entangled state "equal" to this one
        if ent is not None and len(ent) == len(self):
            self._qbits[0]._ent._prepare(self._qbits, lambda x,i: x[i]._pos)
            return self._qbits[0]._ent._state
        state = ket(1)
        for q in self._qbits:
            state = state @ q._get_state()
        return state

    def __matmul__(self, q):
        if isinstance(q, int):
            raise IllegalOperationError("Entangled quantum states cannot be 'duplicated'")
        return Qbits(self._qbits + [q])

    def apply(self, operator, pos=None):
        if pos is not None or not isinstance(operator, (Op, MeasureOp)):
            operator = complete_Op(self._qbits, operator, pos)

        if operator._is_sep():
            i=0
            for op in operator._pieces:
                if len(op) == len(self._qbits[i]):
                    self._qbits[i].apply(op)
                    i+=1
                else:
                    l=len(op)
                    ent=[]
                    while l > 0:
                        ent.append(self._qbits[i])
                        l-=len(self._qbits[i])
                        i+=1
                    qent = gen_ent(ent)
                    qent._apply_qs(op, ent)

        else:
            qent = gen_ent(self._qbits)
            qent._apply_qs(operator, self._qbits)

    def measure(self, basis = None):
        qent = get_ent(self._qbits)
        if qent is None:
            self._state = self._get_state()
            i, state, basis = super().measure(basis)
            self._state = None

            self.apply(MeasureOp(state, basis, i))
            return basis.ew[i]
        return qent.measure(basis)

    def print_as(self, basis):
        is_ent = all((not is_ent(q) for q in self._qbits))
        state = basis.transform(self._get_state())
        s = ''
        for i in range(len(self)):
            if not equal(state[i], 0):
                s += val2str(state[i]) + '|' + basis.symbols[i] + '> '

        if is_ent:
            s = s.replace('+', '±').replace('-', '±')

        return s[:-1]

    @unreal
    def print_probs(self, basis=None):
        b = basis if basis is not None else self._basis
        if b is None:
            T = Matrix((1))
            for q in self._qbits:
                T @= q._basis
            state = rounded_vector(T * self._get_state())
            symb = Vdigit([q._basis.symbols for q in self._qbits])

        else:
            state = rounded_vector(b.transform(self._get_state()))
            symb = b.symbols

        print(format_probs(state, symb))

    @unreal
    def __str__(self):
        if self._basis is None:
            T = Matrix((1))
            ent = False
            for q in self._qbits:
                T @= q.basis
                if is_ent(q): ent = True

            state = rounded_vector(T*self._get_state())
            s = ""
            st = Vdigit([q._basis.symbols for q in self._qbits])

            for e in state:
                val = next(st)
                if e != 0: s += val2str(e) + '|' + val + '> '

            if ent and get_ent(self._qbits) is None:
                s = s.replace('+', '±').replace('-', '±')

            return s[:-1]

        return self.print_as(self._basis)

# ↑↑↑↑↑↑↑↑↑↑↑↑ Qbits class ↑↑↑↑↑↑↑↑↑↑↑↑ #


# ↓↓↓↓↓↓↓↓↓↓↓↓ Qent class ↓↓↓↓↓↓↓↓↓↓↓↓ #

class Qent(Qbits, _Qstate):

    # A particular class that manage the entangled states

    def __init__(self, qbits):
        try:
            qs_list = []
            state = ket(1)

            for q in qbits:
                if q not in qs_list:
                    if is_ent(q) and q._ent is not self:
                        state @= q._ent._state
                        for qent in q._ent._qbits:
                            qent._entangle(self, len(qs_list))
                            qs_list.append(qent)
                    else:
                        state @= q._state
                        q._entangle(self, len(qs_list))
                        qs_list.append(q)

            _Qstate.__init__(self, state, None, qs_list, len(state))

        except Exception as e:
            raise InitializationError('Error to initialize the Qbits entangled', e) from e

    def _get_state(self):
        return self._state

    def apply(self, operator, pos=None):
        if pos is not None:
            index = self._prepare(pos)
            operator = complete_Op(self._qbits, operator, index)

        _Qstate.apply(self, operator)

    def apply2all(self, op):
        self.apply(op @ self.spaces)

    def _apply_qs(self, operator, qs):
        index = self._prepare(qs, get=lambda x, i: x[i]._pos)
        operator = complete_Op(self._qbits, operator, index)
        self.apply(operator)

    def measure(self, basis=None):
        i, _, basis = _Qstate.measure(self, basis)
        self._state = basis[i]
        return basis.ew[i]

    def _prepare(self, pos, get = lambda x,i: x[i]):
        if isinstance(pos, (list,tuple)):
            index = i = min((get(pos, i) for i in range(len(pos))))
            for p in range(len(pos)):
                j = get(pos, p)
                while i > j:
                    j = get(pos, j)
                self._swap(i, j)
                i+=1
            return index

        return pos

    def _swap(self, i, j):
        if i == j: return
        lens = [len(b) for b in self._qbits]
        b = min(len(self[i]), len(self[j]))
        state = Vdigit(lens)
        for x in range(b):
            for y in range(x+1, b):
                while state.has_next:
                    state.ds[i] = x
                    state.ds[j] = y
                    ind_i = state.extract()
                    state.ds[i] = y
                    state.ds[j] = x
                    ind_j = state.extract()

                    # swap!
                    tmp = self._state[ind_i]
                    self._state[ind_i] = self._state[ind_j]
                    self._state[ind_j] = tmp

                    state._incr_ent(i, j, len(state.ds)-1)

        tmp = self._qbits[i]
        self._qbits[i] = self._qbits[j]
        self._qbits[j] = tmp
        # update
        self[i]._entangle(self, i)
        self[j]._entangle(self, j)

    def _calc_p_state(self, pos):
        # try to build the pos state using the probability to measure it
        state=[]
        for s in range(len(self._qbits[pos])):
            state.append(self._get_module(pos, s))
        return ket(*state)

    def _get_module(self, pos, state):
        # get the value of module of a certain state in a certain position
        d = 1
        for i in range(pos+1): d *= len(self._qbits[i])
        dim = len(self._state)//d
        step = dim*len(self._qbits[pos])
        val = 0
        sign = 0
        for i in range(state*dim, len(self._state), step):
            for j in range(dim):
                sign += self._state[i+j]
                val += mod_square(self._state[i+j])

        return math.sqrt(val)

# ↑↑↑↑↑↑↑↑↑↑↑↑ Qent class ↑↑↑↑↑↑↑↑↑↑↑↑ #
