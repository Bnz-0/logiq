from random import random, sample

from .Qerrors import DimensionError, GenericLogiqError, InitializationError
from .qtils import equal, is_scalar, math, mod_square, np


#### Qmath.py
#
# This file contains 2 classes that wrap the numpy.matrix class: the vector and the matrix class,
# it contains also other mathematical functions used in the classes.
# (In the future, these 2 classes could be refactored and became the subclasses of a superclass "linobj")
#
####


# ↓↓↓↓↓↓↓↓↓↓↓↓ Add the npm() method to the np.matrix class ↓↓↓↓↓↓↓↓↓↓↓↓ #
def _npm(self):
    return self

setattr(np.matrix, 'npm', _npm)
# ↑↑↑↑↑↑↑↑↑↑↑↑ Add the npm() method to the np.matrix class ↑↑↑↑↑↑↑↑↑↑↑↑ #


# ↓↓↓↓↓↓↓↓↓↓↓↓ Kronecker product functions ↓↓↓↓↓↓↓↓↓↓↓↓ #

def kron(m1, m2):
    "Return the Kronecker product between `m1` and `m2`  \n(`m1` and `m2` must be `vector` or `matrix`, otherwise they can be a `np.matrix`)"
    try:
        npm1 = m1.npm()
        npm2 = m2.npm()
    except AttributeError as e:
        raise TypeError("The types of the inputs must be vector or matrix") from e
    return select_type(np.kron(npm1, npm2))


def nkron(m, n):
    "Apply the Kronecker product `n` times: `(m @ m @ ... @ m)` for `n` times"
    if n < 1: return Matrix((1)) # [[1]] @ x = x
    out = m
    for _ in range(n-1):
        out = kron(out, m)
    return out

# ↑↑↑↑↑↑↑↑↑↑↑↑ Kronecker product functions ↑↑↑↑↑↑↑↑↑↑↑↑ #


# ↓↓↓↓↓↓↓↓↓↓↓↓ Other useful functions ↓↓↓↓↓↓↓↓↓↓↓↓ #

def select_type(item):
    if len(item.shape) == 0:
        return item
    elif item.shape == (1,1):
        return item[0,0] # scalar
    elif min(item.shape) == 1:
        return Vector(item, copy=False) # vector
    else:
        return Matrix(item, copy=False) # matrix


class npmath:
    # Methods to use, safely, standard operators between scalars, vector, matrix and numpy.matrix

    @staticmethod
    def safe_npm(x):
        try:
            return x if is_scalar(x) else x.npm()
        except:
            raise TypeError("'"+str(type(x))+"' is an invalid type")

    @staticmethod
    def _fun(x, y, fun):
        return select_type(fun(npmath.safe_npm(x), npmath.safe_npm(y)))

    @staticmethod
    def add(x, y):
        return npmath._fun(x, y, lambda x, y: x+y)

    @staticmethod
    def sub(x, y):
        return npmath._fun(x, y, lambda x, y: x-y)

    @staticmethod
    def mul(x, y):
        return npmath._fun(x, y, lambda x, y: x*y)

    @staticmethod
    def div(x, y):
        return npmath._fun(x, y, lambda x, y: x/y)

    @staticmethod
    def equal(x, y):
        if is_scalar(x) or is_scalar(y):
            return x == y
        try:
            x = x.npm()
            y = y.npm()
        except:
            return False
        return npmath.np_equal(x, y)

    @staticmethod
    def np_equal(x, y):
        if x.shape != y.shape:
            return False
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                if x[i,j] != y[i,j]:
                    return False
        return True

# ↑↑↑↑↑↑↑↑↑↑↑↑ Other useful functions ↑↑↑↑↑↑↑↑↑↑↑↑ #


# ↓↓↓↓↓↓↓↓↓↓↓↓ vector class ↓↓↓↓↓↓↓↓↓↓↓↓ #

class Vector:
    """
    + `v`: describes the value of this vector, it can be every kind of object allowed from `numpy.matrix` constructor, or another vector (in that case it will be copied)
    + `normalize` (optional): if True this vector will normalize (i.e. his norm became 1)
    + `values2round` (optional): a list where the value inside will be rounded for example if `values2round=[1]` and the vector is `[1.00001, 12.00001]` it becomes `[1, 12.00001]` (this according to the global precision)
    """

    def __init__(self, v, normalize = False, values2round = None, copy = True):
        try:
            if isinstance(v, (np.matrix, Vector, Matrix)):
                self.vect = v.npm() if not copy else np.matrix(v.npm())
            else:
                self.vect = np.matrix(v, complex)

            if min(self.vect.shape) != 1:
                raise DimensionError('A vector must be have only one dimension')

        except Exception as e:
            raise InitializationError('Error to initialize the vector', e) from e

        if normalize: self.normalize()
        if values2round is not None: self.round_error(values2round)

    def npm(self):
        """
        Return the reference of the `numpy.matrix` associated
        Watch out: it return the reference of the vector (for performance reason), so if you want to modify it first copy it.
        """
        return self.vect

    def is_row(self):
        "`True` if the vector is a 'row vector'"
        return self.vect.shape[0] == 1

    def is_col(self):
        "`True` if the vector is a 'column vector'"
        return self.vect.shape[1] == 1

    def round_error(self, values):
        """Round the error (using `equals()`) of values in `values`.\n
        For example if `values = (0,1)` and this vector is `|5, 0.9999998>` it may became `|5, 1>`"""
        for i in range(len(self)):
            for v in values:
                if equal(self[i], v):
                    self[i] = v
                elif equal(self[i].real, 0) or equal(self[i].imag, 0):
                    self[i] = complex(0 if equal(self[i].real, 0) else self[i].real, 0 if equal(self[i].imag, 0) else self[i].imag)

    def norm(self):
        "Return the norm of this vector, defined as `sqrt(<v|v>)`"
        return math.sqrt(sum(mod_square(x) for x in self))

    def normalize(self):
        "Normalize this vector (i.e. his norm became 1)"
        d = self.norm()
        if d==0: raise DimensionError("The null vector isn't normalizable")
        for i in range(len(self)):
            self[i] /= d
        return self

    def __pos__(self):
        return self

    def __neg__(self):
        return Vector(self.vect.__neg__(), copy=False)

    def __len__(self):
        return max(self.vect.shape)

    def __getitem__(self, i):
        try:
            item = self.vect[0, i] if self.is_row() else self.vect[i, 0]
            return Vector(item) if isinstance(i, slice) else item
        except IndexError as e:
            raise IndexError("Index out of bound (i="+str(i)+" and vector's length = "+str(len(self))+")") from e

    def __setitem__(self, i, value):
        try:
            if self.is_row():
                self.vect[0, i] = value
            else:
                self.vect[i, 0] = value
        except IndexError as e:
            raise IndexError("Index out of bound (i="+str(i)+" and vector's length = "+str(len(self))+")") from e

    def _vector_op(self, v, fun):
        # BUG: numpy.complex128 + v = np.array
        if isinstance(v, Vector):
            if v.vect.shape != self.vect.shape: raise DimensionError("The vector's lengths must be equal")
            return Vector(fun(self.vect, v.vect), copy=False)
        if is_scalar(v):
            return Vector(fun(self.vect, v), copy=False)
        raise TypeError('Other member must be a vector or a scalar')

    def __add__(self, v):
        try:
            return self._vector_op(v, lambda x,y: x+y)
        except Exception as e:
            raise GenericLogiqError('Unable to sum', e) from e

    def __radd__(self, v):
        try:
            return self._vector_op(v, lambda x,y: y+x)
        except Exception as e:
            raise GenericLogiqError('Unable to sum', e) from e

    def __sub__(self, v):
        try:
            return self._vector_op(v, lambda x,y: x-y)
        except Exception as e:
            raise GenericLogiqError('Unable to subtract', e) from e

    def __rsub__(self, v):
        try:
            return self._vector_op(v, lambda x,y: y-x)
        except Exception as e:
            raise GenericLogiqError('Unable to subtract', e) from e

    def __mul__(self, v):
        npv = npmath.safe_npm(v)
        return select_type(self.vect * npv)

    def __rmul__(self, v):
        npv = npmath.safe_npm(v)
        return select_type(npv * self.vect)

    def __truediv__(self, v):
        npv = npmath.safe_npm(v)
        return select_type(self.vect / npv)

    def __rtruediv__(self, v):
        npv = npmath.safe_npm(v)
        return select_type(npv / self.vect)

    def __matmul__(self, v):
        if isinstance(v, int):
            return nkron(self, v)
        return kron(self, v)

    def __eq__(self, v):
        if isinstance(v, (Vector, Matrix, np.matrix)):
            return npmath.np_equal(self.vect, v.npm())
        return False

    def __invert__(self):  # conjugate transpose
        return Vector(self.vect.H)

    def transpose(self):
        "Transform this vector into its transpose"
        self.vect = self.vect.T

    def conj(self):
        "Transform this vector into its conjugate transpose"
        self.vect = self.vect.H

    def __getattr__(self, name):
        if name in 'tT': return Vector(self.vect.T)
        elif name in 'hH': return ~self
        elif name == 'shape': return self.vect.shape
        raise AttributeError("'vector' object has no attribute '"+name+"'")

    def __hash__(self):
        return hash(str(self.vect))

    def __str__(self):
        braket = "<|>"
        i = self.is_col()
        return braket[i] + "; ".join(str(self[i]) for i in range(len(self))) + braket[i+1]

    def __repr__(self):
        return str(self)

    @staticmethod
    def random(dim):
        "Return a random complex vector long `dim`"
        rsign = lambda: sample((1,-1), 1)[0]
        return Vector([complex(rsign()*random(), rsign()*random()) for _ in range(dim)])

# ↑↑↑↑↑↑↑↑↑↑↑↑ vector class ↑↑↑↑↑↑↑↑↑↑↑↑ #


# ↓↓↓↓↓↓↓↓↓↓↓↓ Useful vector constructors ↓↓↓↓↓↓↓↓↓↓↓↓ #

def rounded_vector(v, values2round=(0, 1, -1)):
    "Generate a vector with the values in `values2round` rounded"
    return Vector(v, values2round=values2round)


class ket(Vector):
    """A shortcut to create a column vector

    Usage: `ket(x0, x1, ..., xn)`"""

    def __init__(self, *values):
        if len(values) == 1: values = values[0]
        super().__init__(values)
        if self.is_row():
            if isinstance(values, Vector):
                self.conj() # transformation from ket to bra
            else:
                self.transpose() # no transformation, only creating a vector column


class bra(Vector):
    """A shortcut to create a row vector

    Usage: `bra(x0, x1, ..., xn)`"""

    def __init__(self, *values):
        if len(values) == 1: values = values[0]
        super().__init__(values)
        if self.is_col():
            if isinstance(values, Vector):
                self.conj() # transformation from ket to bra
            else:
                self.transpose() # no transformation, only creating a vector column

# ↑↑↑↑↑↑↑↑↑↑↑↑ Useful vector constructors ↑↑↑↑↑↑↑↑↑↑↑↑ #


# ↓↓↓↓↓↓↓↓↓↓↓↓ matrix class ↓↓↓↓↓↓↓↓↓↓↓↓ #

class Matrix:
    """
    + `M`: describes the value of this vector, it can be every kind of object allowed to `numpy.matrix` constructor, or another matrix (in that case it will be copied)
    + `values2round` (optional): a list where the value inside will be rounded, for example if `values2round=[1]` and the matrix is `[[1.00001, 12.00001],[0.3,1.12]]` it becomes `[[1, 12.00001],[0.3,1.12]]` (this according to the global precision)
    """

    def __init__(self, M, values2round = None, copy = True):
        try:
            if isinstance(M, (np.matrix, Matrix, Vector)):
                self.mtx = M.npm() if not copy else np.matrix(M.npm())
            else:
                try:
                    if isinstance(M[0], ket):
                        M = [[M[j][i] for j in range(len(M))] for i in range(len(M[0]))]
                except: pass
                finally:
                    self.mtx = np.matrix(M, complex)
        except Exception as e:
            raise InitializationError('Error to initialize the matrix', e) from e

        if values2round is not None:
            self.round_error(values2round)

    def round_error(self, values):
        """Round the error (using `equals()`) of values in `values`.
        For example if `values = (0,1)` and this matrix is `[[5, 0.9999998],[-0.000001, 0.1]]` it may became `[[5, 1],[0, 0.1]]`"""
        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                for v in values:
                    if equal(self[i, j], v):
                        self[i, j] = v

    def is_unitary(self):
        "`True` if this matrix is unitary"
        return self.mtx.shape[0] % self.mtx.shape[1] == 0 and np.allclose(np.eye(self.mtx.shape[0]), self.mtx.H * self.mtx)

    def is_orthonormal(self):
        "`True` if this matrix is orthonormal"
        _, r = np.linalg.qr(self.mtx)
        for i in range(len(self)):
            for j in range(i, len(self)):
                if i==j: continue
                elif not equal(r[i,j], 0):
                    return False
        return True

    def nomr(self):
        "return the norm of this matrix"
        return np.linalg.norm(self.mtx)

    def det(self):
        "Return the determinant of this matrix"
        return np.linalg.det(self.mtx)

    def npm(self):
        """
        Return the reference of the `numpy.matrix` associated.
        Watch out: it return the reference of the matrix (for performance reason), so if you want to modify it first copy it.
        """
        return self.mtx

    def __pos__(self):
        return self

    def __neg__(self):
        return Matrix(self.mtx.__neg__(), copy=False)

    def __getitem__(self, i):
        try:
            return select_type(self.mtx[i])
        except IndexError as e:
            raise IndexError("Index out of bound (request="+str(i)+", shape="+str(self.shape)+")") from e

    def __setitem__(self, i, value):
        self.mtx[i] = value

    def _matrix_op(self, M, fun):
        # BUG: numpy.complex128 + M = np.array
        if isinstance(M, Matrix):
            if M.shape != self.shape:
                raise DimensionError("The matrices' shapes must be equal")
            return Matrix(fun(self.mtx, M.mtx), copy=False)
        elif is_scalar(M):
            return Matrix(fun(self.mtx, M), copy=False)
        else:
            raise TypeError('Other member must be a matrix or a scalar')

    def __add__(self, M):
        try:
            return self._matrix_op(M, lambda x,y: x+y)
        except Exception as e:
            raise GenericLogiqError('Unable to add', e) from e

    def __radd__(self, M):
        try:
            return self._matrix_op(M, lambda x,y: y+x)
        except Exception as e:
            raise GenericLogiqError('Unable to add', e) from e

    def __sub__(self, M):
        try:
            return self._matrix_op(M, lambda x,y: x-y)
        except Exception as e:
            raise GenericLogiqError('Unable to subtract', e) from e

    def __rsub__(self, M):
        try:
            return self._matrix_op(M, lambda x,y: y-x)
        except Exception as e:
            raise GenericLogiqError('Unable to subtract', e) from e

    def __mul__(self, M):
        npm = npmath.safe_npm(M)
        return select_type(self.mtx * npm)

    def __rmul__(self, M):
        npm = npmath.safe_npm(M)
        return select_type(npm * self.mtx)

    def __truediv__(self, M):
        npm = npmath.safe_npm(M)
        return select_type(self.mtx / npm)

    def __rtruediv__(self, M):
        npm = npmath.safe_npm(M)
        return select_type(npm / self.mtx)

    def __matmul__(self, M):
        if isinstance(M, int): return nkron(self, M)
        return kron(self, M)

    def __eq__(self, M):
        if isinstance(M, (Vector, Matrix, np.matrix)):
            return npmath.np_equal(self.npm(), M.npm())
        else:
            return False

    def __invert__(self):
        return Matrix(self.mtx.H)

    def transpose(self):
        "Transform this matrix into its transpose"
        self.mtx = self.mtx.T

    def conj(self):
        "Transform this matrix into its conjugate transpose"
        self.mtx = self.mtx.H

    def __getattr__(self, name):
        if name in 'tT': return Matrix(self.mtx.T)
        elif name in 'hH': return ~self
        elif name == 'shape': return self.mtx.shape
        raise AttributeError("'matrix' object has no attribute '"+name+"'")

    def __hash__(self):
        return hash(str(self.mtx)) # fast hash

    def __str__(self):
        return str(self.mtx)

    def __repr__(self):
        return str(self)

    @staticmethod
    def Id(n):
        "Return an `n X n` identity matrix"
        return Matrix(np.identity(n), copy=False)

    @staticmethod
    def filled(shape, val=0):
        """Creates a matrix `shape[0] X shape[1]` filled with the value in `val`
        NB: if shape is a natural number, the matrix will be a square matrix `shape X shape`"""
        if isinstance(shape, int):
            shape = (shape, shape)
        if val == 0:
            return Matrix(np.zeros(shape, complex), copy=False)
        elif val == 1:
            return Matrix(np.ones(shape, complex), copy=False)
        else:
            M = Matrix(np.empty(shape, complex), copy=False)
            for i in range(M.shape[0]):
                for j in range(M.shape[0]):
                    M[i,j] = val
            return M

    @staticmethod
    def random(shape):
        "Return a random complex matrix according to the given shape"
        if isinstance(shape, int): shape = (shape, shape)

        rsign = lambda: sample((1,-1), 1)[0]
        return Matrix([[complex(rsign()*random(), rsign()*random()) for _ in range(shape[0])] for _ in range(shape[1])])

    @staticmethod
    def rand_unitary(dim):
        # (from scipy)
        random_state = np.random
        z = 1/math.sqrt(2)*(random_state.normal(size=(dim, dim)) + 1j*random_state.normal(size=(dim, dim)))
        q, r = np.linalg.qr(z)
        d = r.diagonal()
        q *= d/abs(d)
        return Matrix(q, copy=False)

    @staticmethod
    def rand_orthonormal(dim):
        # (from scipy)
        random_state = np.random
        H = np.eye(dim)
        D = np.ones((dim,))
        for n in range(1, dim):
            x = random_state.normal(size=(dim-n+1,))
            D[n-1] = np.sign(x[0])
            x[0] -= D[n-1]*np.sqrt((x*x).sum())
            # Householder transformation
            Hx = (np.eye(dim-n+1) - 2.*np.outer(x, x)/(x*x).sum())
            mat = np.eye(dim)
            mat[n-1:, n-1:] = Hx
            H = np.dot(H, mat)
            # Fix the last sign such that the determinant is 1
        D[-1] = (-1)**(1-(dim % 2))*D.prod()
        # Equivalent to np.dot(np.diag(D), H) but faster, apparently
        H = (D*H.T).T
        return Matrix(H, copy=False)


# ↑↑↑↑↑↑↑↑↑↑↑↑ matrix class ↑↑↑↑↑↑↑↑↑↑↑↑ #
