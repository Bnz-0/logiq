from .Qerrors import InitializationError, NotAllowError
from .Qmath import kron, math, matrix, np, roundedVector, vector
from .qtils import STD_SYMBOLS


#### Basis.py
#
# This file contains the Basis class and the creation of most used bases:
# - stdbasis: {|0>, |1>}
# - hadamard: {|+>, |->}
# - bell: {|Φ+>, |Ψ+>, |Ψ->, |Φ->}
#
####


# ↓↓↓↓↓↓↓↓↓↓↓↓ Basis classes ↓↓↓↓↓↓↓↓↓↓↓↓ #

class Basis(matrix):
    """
    A Basis serves to measuring and representing a quantum state

    + `basis`: a matrix that describe it or a list of eigenstate
    + `symbols` (optional): the symbols to represent the autosate of this new basis
    """

    def __init__(self, basis, symbols = None, no_cpy = False):
        try:
            if isinstance(basis, Basis):
                super().__init__(basis, no_cpy=True)
                self.symbols = basis.symbols
                self.ew = basis.ew
                return

            super().__init__(basis, no_cpy=no_cpy)

            if (not self.isOrthonormal()):
                raise ValueError("The vectors of the basis must be 'orthonormal' with each other")
            
            self.ew = tuple(np.linalg.eig(self.mtx)[0])

            if symbols != None:
                self.symbols = symbols
                if len(symbols) != len(self):
                    raise ValueError('Length of symbols and basis must be the same')
            elif len(self) <= len(STD_SYMBOLS):
                self.symbols = STD_SYMBOLS[:len(self)]
            else:
                self.symbols = dynSymb(len(self))

        except Exception as e:
            raise InitializationError('Error to initialize Basis', e)
    

    def eigenstate(self, i):
        "returns the i-th eigenstate"       
        if isinstance(i, slice):
            return tuple(self.eigenstate(index) for index in range(*i.indices(len(self))))

        return super().__getitem__((slice(None), i))
    

    def __getitem__(self, i):
        return self.eigenstate(i)


    def __setitem__(self, i, value):
        raise NotAllowError('Basis is an immutable object')
    

    def transform(self, vect):
        """
        Transform the vector `vect` into the one "saw from this basis"

        (i.e. apply the matrix that describe this basis to `vect`)
        """
        return roundedVector(self.mtx*vect if vect.isCol() else vect*self.T)


    def __or__(self, q):
        return q.measure(self)


    def __xor__(self, qs):
        return [q.measure(self) for q in qs]


    def __matmul__(self, other):
        if isinstance(other, Basis):
            return Basis(kron(self, other), no_cpy=True)
        return super().__matmul__(other)

    
    def __invert__(self):
         return self
    

    def conj(self):
        raise NotAllowError('Basis is an immutable object')
    

    def transpose(self):
        raise NotAllowError('Basis is an immutable object')


    def __getattr__(self, name):
        if name in 'tT': return self
        elif name in 'hH': return self
        elif name == 'shape': return self.mtx.shape
        raise AttributeError("'Basis' object has no attribute '"+name+"'")

    
    def __len__(self):
        return self.shape[0]


    def __str__(self):
        return "".join( (f"|{self.symbols[i]}>: {str(self[i])}\n" for i in range(len(self))) )


    def __repr__(self):
        return str(self)


    @staticmethod
    def random(dim, symbols=None):
        "Return a random basis"
        return Basis(matrix.rand_orthonormal(dim), symbols=symbols, no_cpy=True)



class CanonBasis(Basis):
    """
    CanonBasis is a canonical basis (represented by an identical matrix)

    + `dim`: the dimension of this basis
    + `symbols` (optional): the symbols to represent the autosates of this new basis
    """

    def __init__(self, dim, symbols = None):
        if dim < 2:
            raise ValueError('Minimum length allow for a Basis is 2')

        super().__init__(np.identity(dim), symbols, no_cpy=True)
    
    
    def transform(self, vect: vector):
        return vect



class dynSymb:
    # simple class to provide illimitate symbols for big bases, without allocate useless memory

    def __init__(self, length):
        self.length = length


    def __getitem__(self, i):
        if abs(i)<self.length:
            return str(i%self.length)+' '
        raise IndexError() #TODO
    

    def __len__(self):
        return self.length

    def __str__(self):
        return "".join(self)


# ↑↑↑↑↑↑↑↑↑↑↑↑ Basis classes ↑↑↑↑↑↑↑↑↑↑↑↑ #



# ↓↓↓↓↓↓↓↓↓↓↓↓ Creation of most used bases ↓↓↓↓↓↓↓↓↓↓↓↓ #

rd = 1/math.sqrt(2) #reciprocal diagonal ;)

stdbasis = CanonBasis(2)

hadamard = Basis(
    ((rd,rd),
    (rd,-rd)),
    '+-')

bell = Basis(
    ((rd,0,0,rd),
    (0,rd,rd,0),
    (0,rd,-rd,0),
    (rd,0,0,-rd)),
    ['Φ+','Ψ+','Ψ-','Φ-'])

# ↑↑↑↑↑↑↑↑↑↑↑↑ Creation of most used bases ↑↑↑↑↑↑↑↑↑↑↑↑ #
