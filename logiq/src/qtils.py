import math
import re
from collections import defaultdict
from numbers import Complex

import numpy as np


#### qtils.py
#
# This file contains some utility functions, in particular:
# - to manage the precision and provide functions that kepp in mind it
# - to formatting the qtrings of qubit representation
# - the Vdigit class
# - other
#
####


# standard symbols used for representations
STD_SYMBOLS = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'


def prod(l):
    p=1
    for e in l: p *= e
    return p


def find(l, e):
    for i in range(len(l)):
        if e == l[i]: return i
    return None


# ↓↓↓↓↓↓↓↓↓↓↓↓ Precision ↓↓↓↓↓↓↓↓↓↓↓↓ #

current_prec = DEFAULT_PRECISION = np.finfo(np.complex128).precision -1


def set_precision(p = None):
    "Set the precision used by the function `equal()`"
    global current_prec
    if p is None: current_prec = DEFAULT_PRECISION
    elif p > 0:
        current_prec = p
    else:
        raise ValueError('The value of precision must be an int greater than 0')


def what_precision():
    "Returns the current precision"
    return current_prec

# ↑↑↑↑↑↑↑↑↑↑↑↑ Precision ↑↑↑↑↑↑↑↑↑↑↑↑ #



# ↓↓↓↓↓↓↓↓↓↓↓↓ Necessary mathematic functions ↓↓↓↓↓↓↓↓↓↓↓↓ #

def mod_square(c):
    "Returns the module square of a complex number, which is defined as c•c*, where c* is the complex conjugate of `c`"
    return (c * c.conjugate()).real


def equal(n1, n2):
    "Check if 2 complex numbers are equal, keeping in mind possible error of approximation."
    if isinstance(n1, np.inexact):
        return np.around(n1, current_prec) == np.around(n2, current_prec)
    else:
        return n1 == n2


def equals(v1, v2):
    "Check if 2 linear objects are equal, keeping in mind possible error of approximation"
    try: npm1 = v1.npm() ; npm2 = v2.npm()
    except AttributeError: return v1 == v2
    if npm1.shape != npm2.shape: return False
    for i in range(npm1.shape[0]):
        for j in range(npm1.shape[1]):
            if not equal(npm1[i,j], npm2[i,j]): return False
    return True


def isScalar(var):
    "True if `var` is a scalar (i.e. a value in the complexes)"
    return isinstance(var, Complex)

# ↑↑↑↑↑↑↑↑↑↑↑↑ Necessary mathematic functions ↑↑↑↑↑↑↑↑↑↑↑↑ #



# ↓↓↓↓↓↓↓↓↓↓↓↓ To-string utilities ↓↓↓↓↓↓↓↓↓↓↓↓ #

num_format = "{:+.5}"


def set_n_digits(n):
    """
    Sets the number of digits "after the dot" you want to print.  
    For example if you call `set_n_digits(3)` and then `print(my_qubit)`  
    the result will be something like `'+0.707|0> +0.707|1>'`
    """
    if not isinstance(n, int): raise TypeError('The number of digits must be an int')
    if n<1: raise ValueError('Number of digits too low')
    global num_format
    num_format = "{:+."+str(n)+"}"


def round_float(strx):
    # round a float nubmer according to num_format
    for i in range(len(strx)-1, 0, -1):
        if strx[i] != '0':
            return strx[0:i] if strx[i]=='.' else strx[0:i+1]


def val2str(x):
    # generate a string of a complex value according to the number of digits
    if isinstance(x, complex):
        if equal(x.imag, 0):
            return round_float(num_format.format(x.real))
        elif equal(x.real, 0):
            return round_float(num_format.format(x.imag))+'j'
        else:
            return '+('+round_float(num_format.format(x.real)) + round_float(num_format.format(x.imag))+'j)'
    else:
        return round_float(num_format.format(x))


def formatProbs(state, symb):
    # used to formatting the probabilities of a qubit
    format_perc = lambda x: val2str(mod_square(x)*100).replace('+','')+'%'
    s=''
    for i, b in enumerate(symb):
        s += '|' + b + '>: ' + format_perc(state[i])+'\n'

    return s


# ↑↑↑↑↑↑↑↑↑↑↑↑ To-string utilities ↑↑↑↑↑↑↑↑↑↑↑↑ #



# ↓↓↓↓↓↓↓↓↓↓↓↓ From-string utilities ↓↓↓↓↓↓↓↓↓↓↓↓ #

regex = re.compile(r'([+\-]?)\s*([+\-0-9.j]*)\s*\|([^>]+)>')

def str2states(s):
    # transform a string into a dictionare {state : value}
    states = defaultdict(int)
    for sp in re.finditer(regex, s):
        #NB: sign=1, value=2, state=3
        if len(sp[2])==0:
            states[sp[3]] += -1 if sp[1]=='-' else 1
        else:
            states[sp[3]] += complex(sp[1]+sp[2])

    return states


def states2list(states, symbols):
    # - calssical usage: symbols = Basis.symbols
    # - multi-qubit usage: symbols = Vdigit([all bases])

    return [states.get(symb, 0) for symb in symbols]


# ↑↑↑↑↑↑↑↑↑↑↑↑ From-string utilities ↑↑↑↑↑↑↑↑↑↑↑↑ #



# ↓↓↓↓↓↓↓↓↓↓↓↓ Vdigit classes ↓↓↓↓↓↓↓↓↓↓↓↓ #


class Vdigit: # Variable Base digits
    """
    Transform a number into a representation digit-base:  
    every digit of this representations has his base and it's independent of other digits.

    Usage: 
    + `bases`: a list that describe the length of basis for each digit (if it's a list of int) and the symbols that digit can assume (if it's a list fo strings)
    + `start_value` (optional): the initial value

    Example:  
    + `Vdigit([1,2,3], 5) -> '012'`  
    + `Vdigit(['abc','01234'], 10) -> 'c0'`
    """

    def __init__(self, bases, start_value = 0):
        self.ds = []
        self.bases = bases
        self.hasNext = True

        if isinstance(bases[0], int):
            self.symb = lambda i: STD_SYMBOLS
            self.dim = lambda i: self.bases[i]
        else:
            self.symb = lambda i: self.bases[i]
            self.dim = lambda i: len(self.bases[i])

        self.convert(start_value)


    def convert(self, value):
        if value<0: raise ValueError('The value must be greater or equal to 0')

        self.ds = []
        i = len(self.bases)-1
        while value:
            self.ds.append( value%self.dim(i) )
            value = value//self.dim(i)
            i -= 1

        self.ds.reverse()

        if len(self.ds) < len(self.bases):
            self.ds = [0 for _ in range(len(self.bases)-len(self.ds))] + self.ds
        elif len(self.ds) > len(self.bases):
            raise ValueError('Value too high')


    def extract(self):
        "Extracts the integer value represented by this Vdigit"
        b_prod=1 ; val=0
        for i in range(len(self.ds)-1, -1, -1):
            val += self.ds[i]*b_prod
            b_prod *= self.dim(i)
        return val


    def __pp(self, i):
        if i < 0:
            self.hasNext = False
        
        self.ds[i]+=1

        if self.ds[i] >= self.dim(i):
            self.ds[i] = 0
            self.__pp(i-1)


    def _pp_ent(self, i, j, x):
        # custom method used to swap an entanged state
        if x < 0:
            self.hasNext = False
            return
        elif x in (i,j):
            x -= 1
        else:
            self.ds[x] += 1
            if self.ds[x] >= self.dim(x):
                self.ds[x]=0
                x -= 1
            else:
                return

        self._pp_ent(i,j,x)


    def pp(self):
        """
        The "plus plus" method: do a "+1" to the value of this Vdigit
        
        If it's is out of digits raises a StopIteration exception
        """
        if not self.hasNext:
            raise StopIteration

        self.__pp(len(self.bases)-1)
        

    def __iter__(self):
        return self


    def __next__(self):
        out = str(self)
        self.pp() #++
        return out


    def revStr(self):
        """
        Returns the reversed string of the representation  
        (for example if the representation is '1a', this function returns 'a1')
        """
        return "".join(reversed(str(self)))
    

    def __str__(self):
        return "".join( (self.symb(i)[self.ds[i]] for i in range(len(self.bases))) )

# ↑↑↑↑↑↑↑↑↑↑↑↑ Vdigit classes ↑↑↑↑↑↑↑↑↑↑↑↑ #
