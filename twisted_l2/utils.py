"""Defines miscellaneous utility functions."""

from sage.all import *
from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
from .configs import LogOptions

__all__ = ["Ints", "gap_ZZ_to_int", "ga_adjoint", "tietze_to_syllables", "to_indexed_free_group", "load_hap"]

Ints = AbelianGroupGap([0])
_libgap_fail = libgap.eval("fail")

def load_hap(*, tries = 2):
    """
    Loads the GAP package HAP.
    
    Keyword arguments:
    
        - tries: number of attempts.
          (This is needed because on some installations the package does not load on the first try.)
    """
    for j in range(tries):
        try:
            ret = libgap.load_package("HAP")
        except:
            ret = _libgap_fail
            pass
    return ret

def rational_n(q):
    if LogOptions.PRECISION == -1:
        return q
    if not (q in QQ):
        return q
    return round(float(q), LogOptions.PRECISION)

def gap_ZZ_to_int(n):
    """
    Converts an element of AbelianGroupGap([0]) to a Python integer.
    """
    if n not in Ints:
        raise ValueError("n must be an element of AbelianGroupGap([0])")
    return int(n.exponents()[0])

def ga_adjoint(M):
    """
    Computes the adjoint of a matrix with coefficients in Q[G].
    """
    A = M.base_ring()
    G = A.group()
    if M.base_ring() is not G.algebra(QQ):
        raise ValueError("M must be a matrix over a group ring Q[G]")
        
    return matrix(A, M.ncols(), M.nrows(), lambda i,j: M[j,i].antipode())

def tietze_to_syllables(t):
    s = []
    for i in t:
        sgn = sign(i)
        gen = abs(i)-1
        if s:
            s1 = s[-1]
            if s1[0] == gen:
                s1[1] += sgn
            else:
                s.append([gen,sgn])
        else:
            s.append([gen,sgn])
    
    return [tuple(x) for x in s]

def to_indexed_free_group(G, cc):
    F = Groups().free(IntegerRange(G.ngens()))
    QF = F.algebra(QQ)
    cc_free = []

    change_ring = lambda x: QF.sum_of_terms({F(tietze_to_syllables(g.Tietze())):q for g, q in x}, distinct=True)

    for M in cc:
        cc_free.append(matrix(QF, M.nrows(), M.ncols(), lambda i,j: change_ring(M[i,j])))
    return cc_free
