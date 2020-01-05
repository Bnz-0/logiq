from random import random, sample

from .abs_Qstate import _Qstate
from .Basis import Basis, CanonBasis, stdbasis
from .Qbit import Qbit
from .Qbits import Qbits
from .Qerrors import IncomprehensibleStatusError, InitializationError
from .Qmath import ket, vector
from .qtils import STD_SYMBOLS, isScalar, states2list, str2states


#### Qstate.py
#
# This file contains the Qstate class, the "front" class for Qbit, Qbits, ecc...
# and the qbit function, a shortcut to create quantum states
#
####


# ↓↓↓↓↓↓↓↓↓↓↓↓ Qstate calss ↓↓↓↓↓↓↓↓↓↓↓↓ #


class Qstate(_Qstate):
    """
    Qstate is the quantum state class, it serves to create any kind of quantum states:
    
    + qudits (qubit, qutrit, etc), for example |0>
    + composed qudit states, for example |10010>

    Usage:
    
    + `state`: the vector that describes the quantum state
    + `basis` (optional): the basis of this Qstate, it serves to measure the Qstate and to represent it
    + `normalize` (optional): if True the state will be normalized
    + `transform` (optional): if True the state will be interpreted as a basis representation (so if state=ket(1,0) and basis={|+>,|->} the state will be interpreted as |+> and not |0>)
    """

    def __init__(self, state, basis=None, normalize=False, transform=False):
        
        if isinstance(state, vector):
            self.__class__ = Qbit
            self.__init__(state, basis, normalize, transform)
        
        else:
            self.__class__ = Qbits
            self.__init__(state, basis)


    @staticmethod
    def random(n = 2):
        "Generate a random quantum state of size `n`"
        if isinstance(n, Basis):
            b = n
            n = len(b)
        else:
            b = CanonBasis(n)
        rsign = lambda: sample((1,-1), 1)[0]
        return Qbit(ket([complex(rsign()*random(), rsign()*random()) for _ in range(n)]), b, normalize=True)


    @staticmethod
    def parse(s, basis = None, normalize = None):
        """
        Generate a quantum state parsing it from a string `s`
        
        The string must be formatted according to the classical representation of quantum states (bra-ket notation):  
        "+n1|s1> +n2|s2> ..." where `n` are complex numbers and `s` are the states of the basis `basis`
        
        Usage:
        + `s`: the string to be parsed
        + `basis` (optional): the Basis used in the `s` representation
        + `normalize` (optional): if True normalize the resulting state
        """
        try:
            state = str2states(s)
            if basis is None:
                if all((b in stdbasis.symbols for b in state)):
                    basis = stdbasis
                elif all((b in STD_SYMBOLS for b in state)):
                    basis = CanonBasis(STD_SYMBOLS.find(max(state.keys()))+1)
                elif len(state)>1:
                    basis = CanonBasis(len(state), list(state.keys()))
                else:
                    raise IncomprehensibleStatusError('Fail to parse')

            return Qbit(
                    ket(*states2list(state, basis.symbols)),
                    basis,
                    normalize=normalize,
                    transform=True
                )
            
        except Exception as e:
            raise InitializationError("Error to initialize Qbit from string", e)

# ↑↑↑↑↑↑↑↑↑↑↑↑ Qstate class ↑↑↑↑↑↑↑↑↑↑↑↑ #



def qbit(*something, basis=None, normalize=None, transform=False):
    "A shortcut to generate quantum states (see documentation for a complete description)"

    if len(something)>1:
        if isScalar(something[0]):
            return Qstate(ket(*something), basis, normalize=normalize, transform=transform)
        else:
            return Qstate(something, basis=basis)

    something = something[0] #single element
    if isinstance(something, vector):
        return Qstate(something, basis, normalize=normalize, transform=transform)

    elif isinstance(something, (list, tuple)):
        return Qstate(ket(*something), basis, normalize=normalize, transform=transform)

    elif isinstance(something, str):
        return Qstate.parse(something, basis, normalize=normalize)

    else:
        raise IncomprehensibleStatusError('Failed to understand which quantum state to generate')
