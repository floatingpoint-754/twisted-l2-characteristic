from sage.all import *
from sage.groups.group import is_Group
import itertools
from .configs import FGARankOptions as opt
from .logger import Log, SILENT, INFO, DEBUG
from .utils import rational_n

__all__ = ["finite_group_algebra_rank"]

def try_invert(L, x):
    """
    Tries to compute the inverse of an element x of the group algebra of L over a field.
    Returns the inverse of x if x is recognized as a unit, else None.
    """
    QL = x.parent()
    if L.algebra(QL.base_ring()) is not QL:
        raise ValueError("x must be an element of a group algebra K[L]")
        
    l = list(x)
    if len(l) == 0:
        return None
    if len(l) == 1:
        return l[0][1].inverse_of_unit()*QL(l[0][0].inverse())
    if len(l) == 2:
        a,r = l[0]; b,s = l[1]
        g = a.inverse()*b; q = s/r; m = g.order()
        if 1 == (-q)**m:
            return None
        return 1/(r*(1-(-q)**m)) * sum((-q)**i * QL(g**i) for i in range(m)) * a.inverse()
    return None

def try_reduce(L, M):
    """
    Tries to reduce the matrix M over the group algebra Q[L] using Gauss moves.
    Returns a quadruple (i,j,r,d) such that rk(M) = r + rk(M'[i:,j:])
    and M'[i:,j:]*d has integer coefficients, where M' is M after an invocation of try_reduce.
    Nothing is checked.
    """
    A = L.algebra(QQ)
    if type(M) != sage.matrix.matrix_generic_dense.Matrix_generic_dense or M.base_ring() != A:
        raise ValueError("M must be a matrix over Q[S]")
        
    ii = 0; jj = 0
    r, c = M.dimensions()
    res = 0
    step = False
    # try to find an invertible element somewhere
    while True:
        N = matrix(ZZ, r-ii, c-jj, lambda i,j : len(M[ii+i,jj+j]))
        
        # ===== Get rid of zero rows and cols =====
        
        zerocols = vector([1]*(r-ii)) * N
        zerorows = N * vector([1]*(c-jj))
        freeR = 0; freeC = 0
        for i, el in enumerate(zerorows):
            if el:
                continue
            M.swap_rows(freeR+ii, i+ii)
            N.swap_rows(freeR, i)
            freeR += 1
        
        for j, el in enumerate(zerocols):
            if el:
                continue
            M.swap_columns(freeC+jj, j+jj)
            N.swap_columns(freeC, j)
            freeC += 1
        ii += freeR
        jj += freeC
        
        # ===== try to find some position with a trivial unit, if there's none then a bicyclic unit =====
        
        idx = next(filter(lambda x:N[x] == 1, itertools.product(range(freeR, N.nrows()), range(freeC, N.ncols()))), None)
        if idx is None and opt.BICYCLIC:
            idx = next(filter(lambda x:N[x] == 2, itertools.product(range(freeR, N.nrows()), range(freeC, N.ncols()))), None)
        
        if idx is not None:
            i = idx[0] + ii-freeR; j = idx[1] + jj-freeC
            #print(len(M[i,j]), end=' ')
            x = try_invert(L, M[i,j])
            if x is not None:
                M.swap_rows(i, ii)
                M.swap_columns(j, jj)
                M.set_row_to_multiple_of_row(ii, ii, x)
                for k in range(ii+1,r):
                    M.add_multiple_of_row(k, ii, -M[k, jj])
                res += 1
                ii += 1
                jj += 1
                step = True
            else:
                Log(DEBUG, f"Found no inverse for unit M[{i},{j}]")
        else:
            Log(DEBUG, "Found no easily invertible unit")
            Log(DEBUG, "Lengths of remaining matrix entries:")
            Log(DEBUG, N[freeR:,freeC:])
            Log(DEBUG, "")
        
        if not step:
            break
        step = False
    
    # ===== common denominator =====
    denom = lcm(q.denominator() for i in range(ii, M.nrows()) for j in range(jj, M.ncols()) for g,q in M[i,j])
    
    return (ii, jj, res, denom)

def finite_group_algebra_rank(L, M):
    """
    Computes the von Neumann rank of a matrix over the group ring of a finite group.
    
    Arguments:
    
        - L: finite group.
        - M: matrix over Q[L].
    """
    if not (is_Group(L) or L.is_finite()):
        raise ValueError("L must be a finite group")
    A = L.algebra(QQ)
    if type(M) != sage.matrix.matrix_generic_dense.Matrix_generic_dense or M.base_ring() != A:
        raise ValueError("M must be a dense matrix over Q[L]")
    ring = {"rational": ZZ, "real": RR, "double": RDF}[opt.RING]
    if opt.TRANSPOSE:
        M = M.transpose()
    r = M.nrows()
    c = M.ncols()
    if r == 0 or c == 0:
        return 0
    k = 0
    res_rk = 0
    denom = 1
    N = None
    if opt.TRY_REDUCE:
        ii,jj,res_rk, denom = try_reduce(L, M)
        Log(INFO, "")
        Log(INFO, f"Matrix of size {(r, c)} reduced by {(ii, jj)}.")
        Log(DEBUG, f"Need to add {res_rk} to rank")
        Log(DEBUG, f"Denominator of reduced matrix: {denom}")
        r -= ii; c -= jj
        M = M.submatrix(ii,jj)
    elem = L.list()
    k = len(elem)
    N = matrix(ring, r*k, c*k, sparse=opt.SPARSE)
    for i in range(r):
        for ii in range(k):
            for j in range(c):
                if M[i,j] == A.zero():
                    continue
                for jj in range(k):
                    gg = elem[ii]**(-1)*elem[jj]
                    q = M[i,j][gg]
                    if q != 0:
                        q *= denom
                        N[i*k+ii, j*k+jj] = q

    Log(INFO, f"Dimensions of N: {N.dimensions()}")  
    Log(DEBUG, f"Nonzero entries of N: {len(N.support())} out of {N.nrows() * N.ncols()}")
    res = 0
    if ring is RDF:
        sv = N.singular_values("auto")
        res = len(sv) - sv.count(0)
    elif ring is ZZ:
        res = N.rank("linbox" if opt.SPARSE else "flint")
    else:
        res = N.rank()
    total_rank = res_rk + Integer(res) / k
    Log(INFO, f"Rank: {rational_n(total_rank)} (rounds up to {ceil(total_rank)})")
    
    return ceil(total_rank) if opt.CEILING else total_rank
