from .Basis import Basis, CanonBasis, hadamard
from .Qerrors import DimensionError, InitializationError, NotAllowError
from .Qmath import ket, kron, Matrix, nkron, np, npmath, Vector
from .qtils import Vdigit, find, is_scalar, mod_square, states2list, str2states


### Operator.py
#
# This file contains the Op class and the creation of most used operators
#
####


# ↓↓↓↓↓↓↓↓↓↓↓↓ Op classes ↓↓↓↓↓↓↓↓↓↓↓↓ #

class Op(Matrix):
    """
    An Op (operator) server to modify the state of a quantum state

    + `operator`: a matrix that represent the operator
    """

    def __init__(self, operator, _pieces = None, copy = True):
        try:
            if isinstance(operator, Op):
                super().__init__(operator.mtx, copy=False)
                self._pieces = operator._pieces
            else:
                super().__init__(operator, copy=copy)
                if not self.is_unitary():
                    raise ValueError('Matrix not unitary')

            self._pieces = [self] if _pieces is None else _pieces

        except Exception as e:
            raise InitializationError("Error to initialize Op", e) from e

    def _is_sep(self):
        return len(self._pieces) > 1

    def __setitem__(self, i, value):
        raise NotAllowError('Op is an immutable object')

    def __or__(self, q):
        q.apply(self)

    def __xor__(self, q):
        q.apply2all(self)

    def __mul__(self, other):
        if isinstance(other, Op):
            return Op(self.mtx * other.mtx, copy=False)

        if is_scalar(other) and mod_square(other) == 1:
            return Op(self.mtx * other)

        return super().__mul__(other)

    # def __rmul__(self, other):
    #     return super().__rmul__(other)

    def __matmul__(self, other):
        if isinstance(other, Op):
            return Op(kron(self.mtx, other.mtx), _pieces=self._pieces+other._pieces, copy=False)

        if isinstance(other, int):
            return Op(nkron(self, other), _pieces=self._pieces*other, copy=False)

        return super().__matmul__(other)

    def __eq__(self, other):
        return npmath.equal(self, other)

    def __invert__(self):
        return Op(self.mtx.H)

    def __len__(self):
        return self.mtx.shape[0]

    @staticmethod
    def build(rules, basis = None):
        """
        Another way to create an operator is to describe it's behavior with the autostate of a basis

        + `rules`: a dictionary {state : new_state} where state and new_state can be the vector, the number of autostate of basis or a string that represents the state (according to the basis symbols)
        + `basis` (optional): the basis which represents the rules (if it's the standard basis it's unnecessary)
        """

        op = Matrix.filled(len(rules), 0)

        if basis is None:
            basis = CanonBasis(len(rules))
        elif isinstance(basis, (list,tuple)):
            # TODO: improve with __getitem__ in Vdigit && lazy init of CanonBasis
            symb = [str(v) for v in Vdigit(basis)]
            basis = CanonBasis(len(symb), symb)

        for k, v in rules.items():
            i = find(basis.symbols, k.strip('|> ')) if isinstance(k, str) else k

            if isinstance(v, str):
                state = basis.transform(ket(*states2list(str2states(v), basis.symbols)))
            elif isinstance(v, int):
                state = basis[v] if v>=0 else -basis[-v] #NB: beware of state 0: -0 is 0
            else:
                state = v

            op += state * ~basis[i]

        return Op(op.npm(), copy=False)

    @staticmethod
    def neutral():
        "Returns an operator `N` such that `N @ x = x`"
        return Op([1], _pieces=[])
        #[[1]] @ x = x

    @staticmethod
    def Id(n):
        "Return an `n X n` identity operator"
        return Op(np.identity(n), copy=False)

    @staticmethod
    def phase_gate(phase, deg = False):
        """Generates the 'phase gate' according to the giving `phase`
        if `deg` is `True` the phase will be evaluated in degrees and not in radians"""
        if deg: phase = (2*np.pi*phase)/360.0
        return Op([[1, 0], [0, np.e**(1j * phase)]])

    @staticmethod
    def C(U):
        """
        The Controlled-U gate

        + `U`: the operator to apply to the second qubit if the first is "true"
        """
        if len(U) != 2: raise DimensionError('The operator U must be a 2x2 matrix')
        return Op([
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, U[0,0], U[0,1]],
            [0, 0, U[1,0], U[1,1]]
        ])

    @staticmethod
    def random(n):
        "Generates a random `n X n` operator"
        return Op(Matrix.rand_unitary(n), copy=False)


class MeasureOp(Matrix):
    # A special type of operator to permit the creation of a "measurement operator"
    # (because normally an operator can't be non-unitary)

    def __init__(self, state, basis, i):
        super().__init__(((basis[i] * ~basis[i])/state[i]).npm(), copy=False)
        self._pieces = [self]

    def _is_sep(self):
        return False

    def __or__(self, q):
        q.apply(self)

    def __xor__(self, q):
        q.apply2all(self)

    def __len__(self):
        return self.mtx.shape[0]


# ↑↑↑↑↑↑↑↑↑↑↑↑ Op classes ↑↑↑↑↑↑↑↑↑↑↑↑ #


# ↓↓↓↓↓↓↓↓↓↓↓↓ Creation of most used Operators ↓↓↓↓↓↓↓↓↓↓↓↓ #

Op.I = Op( [[1,0],[0,1]] )
Op.X = Op( [[0,1],[1,0]] )
Op.Y = Op( [[0,-1j],[1j,0]] )
Op.Z = Op( [[1,0],[0,-1]] )
Op.H = Op(hadamard.npm(), no_cpy=True)
Op.cnot = Op( [[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]] )
Op.swap = Op( [[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]] )
Op.sqrtSwap = Op( [[1,0,0,0],[0,0.5*(1+1j),0.5*(1-1j),0],[0,0.5*(1-1j),0.5*(1+1j),0],[0,0,0,1]] )

# ↑↑↑↑↑↑↑↑↑↑↑↑ Creation of most used Operators ↑↑↑↑↑↑↑↑↑↑↑↑ #
