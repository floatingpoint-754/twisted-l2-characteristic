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

def Uplift(M):
    """
    (DEPRECATED) Lifts M to the group ring of the free group.
    """
    A = M.base_ring()
    G = A.group()
    R = A.base_ring()
    if M.base_ring() is not G.algebra(R):
        raise ValueError("M must be a matrix over a group ring R[G]")
    F = G.free_group()
    B = F.algebra(R)
    
    up = lambda g: F(g.Tietze())
    upr = lambda x: B.sum_of_terms({up(k): v for k,v in x}, distinct=True)
    return matrix(B, M.nrows(), M.ncols(), lambda i,j: upr(M[i,j]))

def Pushdown(M, G):
    """
    (DEPRECATED) Coerces M from the group ring of the free group.
    (DOES NOT WORK: GAP sometimes attempts to solve the word problem)
    """
    B = M.base_ring()
    F = B.group()
    R = B.base_ring()
    if M.base_ring() is not F.algebra(R):
        raise ValueError("M must be a matrix over a group ring R[F]")
    if G.free_group() is not F:
        raise ValueError("G must be a quotient of F")
    A = G.algebra(R)
    
    down = lambda f: G(f.Tietze())
    def downr(x):
        lst = [v*B.monomial(down(k)) for k,v in x.monomial_coefficients().items()]
        print(lst)
        return sum(lst)
    return matrix(A, M.nrows(), M.ncols(), lambda i,j: downr(M[i,j]))

def Laplacian(M, N):
    """
    Computes the combinatorial Laplacian of two matrices with coefficients in R[F].
    """
    A = M.base_ring()
    F = A.group()
    R = A.base_ring()
    if M.base_ring() is not F.algebra(R) or N.base_ring() is not F.algebra(R):
        raise ValueError("M, N must be matrices over a group ring R[F]")
    if M.ncols() != N.nrows():
        raise ValueError("M, N are not compatible")
    
    Mt = ga_adjoint(M)
    Nt = ga_adjoint(N)
    return Mt*M+N*Nt

def SimplifyAlgebraElement(x):
    """
    (DEPRECATED: Tries to solve the word problem)
    """
    A = x.parent()
    ds = {}
    for k,v in x:
        found = False
        for kk in ds:
            if kk == k:
                ds[kk] += v
                found = True
                break
        if not found:
            ds[k] = v
    
    return A.linear_combination(ds.items())

def AlgebraElementRepr(x):
    coeffs = x.monomial_coefficients()
    if len(coeffs) == 0:
        return "0"
    
    first = True
    res = ""
    for k,v in x.monomial_coefficients().items():
        if not first:
            res += " + " if v > 0 else " - "
            v = v if v > 0 else -v
        res += ((str(v) + "*") if v != 1 else "") + str(k)
        first = False
    return res

def PrintMatrix(M):
    ls = [[AlgebraElementRepr(M[i,j]) for j in range(M.ncols())] for i in range(M.nrows())]
    ln = max([len(s) for s in flatten(ls)], default = 0)
    for lr in ls:
        print("[" + "  ".join([s.rjust(ln) for s in lr]) + "]")    
        
