"""
Defines functions computing the twisted rank of matrices with coefficients in a group ring.
"""

from sage.all import *
from .utils import tietze_to_syllables, Ints, gap_ZZ_to_int, rational_n
from .fga_rank import finite_group_algebra_rank
import functools
from .logger import Log, ProgressBar, SILENT, INFO, DEBUG, PB_START, PB_STOP
from .configs import TwistedRankOptions as opt

__all__ = [
    "finite_quotient",
    "von_neumann_rank",
    #"determinant_degree",
]

@functools.cache
def finite_quotient(G, l):
    exps = [(Primes().unrank(i), e) for i, e in enumerate(l) if e]
    n = len(exps)
    if n == 0:
        G1 = libgap.TrivialGroup()
        epi = libgap.GroupHomomorphismByImagesNC(G, G1, [G1.One()]*G.ngens())
        return (epi, epi.Image())
    if n == 1:
        p,c = exps[0]
        epi = libgap.EpimorphismPGroup(G, p, c)
        return (epi, epi.Image())
    epis = [libgap.EpimorphismPGroup(G, p, c) for p,c in exps]
    imgs = [f.Image() for f in epis]
    img = libgap.DirectProduct(*imgs)
    gens_imgs = [f.MappingGeneratorsImages()[1] for f in epis]
    embed_imgs = [[libgap.Embedding(img, i+1).Image(g) for g in gi] for i, gi in enumerate(gens_imgs)]
    prod_imgs = [libgap.Product([gi[i] for gi in embed_imgs]) for i in range(G.ngens())]
    return (libgap.GroupHomomorphismByImagesNC(G, img, prod_imgs), img)

def von_neumann_rank(G, M, exps, LogObj = lambda: None):
    """
    Computes the von Neumann rank of a matrix with coefficients in the group ring Q[G].
    
    Arguments:
    
        - G: finitely presented group
        - M: matrix over Q[G] acting by right multiplication (given as a matrix over the free group)
        - exps: nilpotency class list for EpimorphismPGroup approximation
        - LogObj: an object provided so that this method can save data as members
    
    Complexity:
    
        The matrix M is enlarged by a factor possibly exponential in p and c.
    """
    
    if G not in Groups():
        raise ValueError("G must be a group")
    
    
    B = M.base_ring()
    Fb = B.group()
    from sage.groups.indexed_free_group import IndexedFreeGroup
    
    k = M.nrows()
    h = M.ncols()

    if k == 0 or h == 0:
        return Integer(0)
    
    # Convert to IndexedFreeGroup, which is faster
    if not isinstance(Fb, IndexedFreeGroup):
        F = Groups().free(IntegerRange(G.ngens()))
        QF = F.algebra(QQ)
        #chfree = lambda x: QF.sum_of_terms({F(tietze_to_syllables(g.Tietze())):q for g, q in x}, distinct=True)
        _M = matrix(QF, k, h)
        for i in range(k):
            for j in range(h):
                _M[i,j] = QF.sum_of_terms(((F(tietze_to_syllables(g.Tietze())), q) for g, q in M[i,j]), distinct=True)
        M = _M
    else:
        F = Fb
    
    gens = {k for i in range(k) for j in range(h) for k,v in M[i,j]}
    gens.discard(F.one())
    fhom, Fin = finite_quotient(G, exps)
    Log(INFO, f"Size of Fin: {Fin.Size()}")
    LogObj.Fin = Fin
    LogObj.Sgens = gens
    
    Log(DEBUG, f"Making finite quotient L... need to compute {len(gens)} images")
    
    genset = set()
    gendict = {F.one(): Fin.One()}
    
    ProgressBar(DEBUG, PB_START)
    for i,gensx in enumerate(gens):
        G_element = G([s*(j+1) for j,s in gensx.to_word_list()])
        imag = fhom.Image(G_element)
        genset.add(imag)
        gendict[gensx] = imag
        ProgressBar(DEBUG, i/len(gens))
    
    ProgressBar(DEBUG, PB_STOP)
    Log(DEBUG, "Calling libgap...")
    
    # construct the smallest group/ring containing the matrix coefficients
    # we use rk_G(N(G) tensor M) = rk_S(M) for S < G
    
    L = Fin.Subgroup(list(genset)) 
    Log(INFO, f"|L| = {L.Size()}")
    LogObj.Lsize = L.Size()
    
    LPhom = L.IsomorphismPermGroup()
    LP = LPhom.Image()
    small = LP.SmallerDegreePermutationRepresentation()
    LL = small.Image()
    LLSage = PermutationGroup(gap_group=LL)
    C = LLSage.algebra(QQ)
    LogObj.L = LL
    
    L_to_LLSage = {}
    def F_to_LLSage(k): # this HALVES the runtime when rank is cheap
        try:
            k_L = gendict[k]   
        except KeyError:
            raise RuntimeError("This should never happen")
        try:
            return L_to_LLSage[k_L]
        except KeyError:
            imag = LLSage(small.Image(LPhom.Image(k_L)))
            L_to_LLSage[k_L] = imag
            return imag
    
    Log(INFO, "Constructing matrix over Q[L]...")
    AF = matrix(C, k, h)
    
    ProgressBar(DEBUG, PB_START)
    for i in range(k):
        for j in range(h):
            AF[i,j] = C.sum_of_terms(((F_to_LLSage(k), v) for k,v in M[i,j]))
            ProgressBar(DEBUG, (i*h + j)/(k*h))
    ProgressBar(DEBUG, PB_STOP)
    
    if opt.PRINT_FIN:
        print(AF)
    if opt.ASK_INPUT and input("Keep a copy of the matrix over Q[L]? "):
        LogObj.finite_alg_matrix = AF
        
    Log(INFO, "Computing rank...")
    rnk = finite_group_algebra_rank(LLSage, AF)
        
    #Log(INFO, f"Rank: {rational_n(rnk)}")
    return rnk

def matrix_expansion(G, phi, lift, M, n, LogObj = lambda: None, uniform_bias = False):
    """
    Expands a matrix with coefficients in Q[G], according to Oki's matrix expansion algorithm (arXiv:1907.04512, 2019).

    Arguments:

        - G: finitely presented group
        - phi: surjective morphism G --> Z (can be given as a list of integers of length G.ngens())
        - lift: an element of the free group of G such that phi(lift) == 1
        - M: matrix over Q[G] (given as a matrix over the free group)
        - n: upper bound for matrix expansion
        - LogObj: an object provided so that this method can save data as members
        - uniform_bias: multiply each row by the same power of lift

    Returns:
        A tuple (Mx, bias), where Mx is the expanded matrix and bias is either
            - the sum of all shifts applied to each row (if not uniform_bias)
            - the shift applied to each row (if uniform_bias)

    Complexity:

        Polynomial in n.
    """

    if G not in Groups():
        raise ValueError("G must be a group")
    if hasattr(phi, "parent") and phi.parent() == Hom(G, Ints):
        phi = [gap_ZZ_to_int(phi(g)) for g in G.gens()]
    elif not isinstance(phi, list) and not isinstance(phi, tuple):
        raise ValueError("phi must be a morphism G --> Z")
    if abs(gcd(phi)) != 1:
        raise ValueError("phi must be surjective")

    x = lift
    B = M.base_ring()
    F = B.group()
    from sage.groups.indexed_free_group import IndexedFreeGroup
    if not isinstance(F, IndexedFreeGroup):
        raise ValueError(f"The base group of M must be an IndexedFreeGroup (got {F})")

    M = copy(M)
    k, h = M.nrows(), M.ncols()

    # First, we multiply the rows on the left by adequate multiples of x,
    # such that all of M's entries have positive valuation

    def val_F(f):
        return sum((s*phi[j] for j,s in f.to_word_list()), 0)

    def valmin(y):
        """Computes the order valuation of a polynomial."""
        return min((val_F(f) for f,q in y), default=0)

    minrows = [min((valmin(M[i,j]) for j in range(h)), default=0) for i in range(k)]
    if not uniform_bias:
        bias = sum(minrows)
        Log(DEBUG, f"Minimum exponent of each row: {minrows}, total: {bias}")
    else:
        bias = min(minrows)
        minrows = [bias] * k
        Log(DEBUG, f"Multiplying each row by u^{-bias}")


    for i in range(k):
        for j in range(h):
            M[i,j] = B.monomial(x**(-minrows[i]))*M[i,j]


    def polynomial_degree(y):
        return max((val_F(f) for f,q in y), default=0)

    max_degree = max((polynomial_degree(entry) for entry in M.list()), default=0)
    if k != h:
        mu = n
        Log(INFO, f"Maximum valuation of entries = {max_degree}, expanding by {mu}")
    else:
        mu = min(n, max_degree*k)
        Log(INFO, f"Maximum valuation of entries = {max_degree}, should expand matrix by {max_degree*k}, expanding by {mu} instead")

    if mu == 0:
        Log(INFO, "Valuation bound: 0 (skipping this...)")
        return matrix(B, 0, 0), bias

    Ax = matrix(B, mu*k, mu*h)

    ProgressBar(DEBUG, PB_START)
    for i in range(k):
        for j in range(h):
            l = [B.monomial(x**s) * M[i,j] for s in range(mu)]
            for s in range(mu):
                for co, q in l[s]:
                    H = val_F(co)
                    if H in range(mu):
                        Ax.add_to_entry(i+k*s, j+h*H, q * B.monomial(co * x**(-H)))

            ProgressBar(DEBUG, (i*h+j)/(k*h))

    ProgressBar(DEBUG, PB_STOP)
    Log(DEBUG, "expanded matrix done")
    return Ax, bias

def determinant_degree(G, phi, lift, M, n, exps, LogObj = lambda: None):
    """
    Computes the degree of the Dieudonné determinant of a self-adjoint square matrix with coefficients in the group ring Q[G].
    
    Arguments:
    
        - G: finitely presented group
        - phi: surjective morphism G --> Z (can be given as a list of integers of length G.ngens())
        - M: self-adjoint square matrix over Q[G] (given as a matrix over the free group)
        - n: upper bound for matrix expansion
        - exps: nilpotency class list for EpimorphismPGroup approximation
        - LogObj: an object provided so that this method can save data as members

    Complexity:
    
        The matrix M is enlarged by a factor linear in n and possibly exponential in p and c.
    """
    
    if G not in Groups():
        raise ValueError("G must be a group")
    if hasattr(phi, "parent") and phi.parent() == Hom(G, Ints):
        phi = [gap_ZZ_to_int(phi(g)) for g in G.gens()]
    elif not isinstance(phi, list) and not isinstance(phi, tuple):
        raise ValueError("phi must be a morphism G --> Z")
    if abs(gcd(phi)) != 1:
        raise ValueError("phi must be surjective")

    x = lift
    LogObj.lift = x
    Log(INFO, f"Lift: {x}")
    
    B = M.base_ring()
    F = B.group()
    from sage.groups.indexed_free_group import IndexedFreeGroup
    if not isinstance(F, IndexedFreeGroup):
        raise ValueError(f"The base group of M must be an IndexedFreeGroup (got {F})")

    k = M.nrows()
    
    Ax, bias = matrix_expansion(G, phi, lift, M, n, LogObj)

    if opt.ASK_INPUT and input("Keep a copy of the expanded matrix? "):
        LogObj.expanded_matrix = copy(Ax)

    rnk = von_neumann_rank(G, Ax, exps, LogObj)
    valn = Ax.nrows() - rnk

    Log(INFO, f"Valuation: {rational_n(valn)}")
    LogObj.valuation = valn
    
    deg = (valn + bias) * -2
    Log(INFO, f"Degree: {rational_n(deg)}")
    return deg

def determinant_degree_asymmetric(G, phi, lift, M, n, exps, LogObj = lambda: None):
    """
    Computes the degree of the Dieudonné determinant of a general square matrix with coefficients in the group ring Q[G].

    Arguments:

        - G: finitely presented group
        - phi: surjective morphism G --> Z (can be given as a list of integers of length G.ngens())
        - M: square matrix over Q[G] (given as a matrix over the free group)
        - n: upper bound for matrix expansion
        - exps: nilpotency class list for EpimorphismPGroup approximation
        - LogObj: an object provided so that this method can save data as members

    Complexity:

        The matrix M is enlarged by a factor linear in n and possibly exponential in p and c.
    """

    if G not in Groups():
        raise ValueError("G must be a group")
    if hasattr(phi, "parent") and phi.parent() == Hom(G, Ints):
        phi = [gap_ZZ_to_int(phi(g)) for g in G.gens()]
    elif not isinstance(phi, list) and not isinstance(phi, tuple):
        raise ValueError("phi must be a morphism G --> Z")
    if abs(gcd(phi)) != 1:
        raise ValueError("phi must be surjective")

    x = lift
    LogObj.lift = x
    Log(INFO, f"Lift: {x}")
    minus_phi = [-z for z in phi]

    B = M.base_ring()
    F = B.group()
    from sage.groups.indexed_free_group import IndexedFreeGroup
    if not isinstance(F, IndexedFreeGroup):
        raise ValueError(f"The base group of M must be an IndexedFreeGroup (got {F})")

    k = M.nrows()

    A0, bias0 = matrix_expansion(G, phi, lift, M, n, LogObj, uniform_bias=True)
    A1, bias1 = matrix_expansion(G, minus_phi, lift**-1, M, n, LogObj, uniform_bias=True)

    if opt.ASK_INPUT and input("Keep a copy of the expanded matrices? "):
        LogObj.expanded_matrix0 = copy(A0)
        LogObj.expanded_matrix1 = copy(A1)

    val0 = A0.nrows() - von_neumann_rank(G, A0, exps, LogObj)
    val1 = A1.nrows() - von_neumann_rank(G, A1, exps, LogObj)

    Log(INFO, f"Valuations: {rational_n(val0)}, {rational_n(val1)}")
    Log(INFO, f"Biases: {rational_n(bias0)}, {rational_n(bias1)}")
    LogObj.valuation = (val0, val1)

    deg = - val0 - bias0 - val1 - bias1
    Log(INFO, f"Degree: {rational_n(deg)}")
    return deg

def dim_torsion_coker(G, phi, lift, M, rk, n, exps, LogObj = lambda: None):
    """
    Given a matrix M with coefficients in the group ring Q[G],
    returns the dimension over D(K) of the torsion D(K)[u, u^-1]-module of coker(M).
    This is equivalent to the sum of the degrees of nonzero entries in the Jacobson normal form of M,
    and generalizes the method `determinant_degree_asymmetric` to singular matrices.

    Note: using `determinant_degree_asymmetric` is still more efficient for nonsingular square matrices,
    as that method allows every row to be multiplied by a different power of the lift.

    The output is not meaningful when k is not the rank of M over D(G).

    Arguments:

        - G: finitely presented group
        - phi: surjective morphism G --> Z (can be given as a list of integers of length G.ngens())
        - M: matrix over Q[G] (given as a matrix over the free group)
        - rk: rank of M over D(G)
        - n: upper bound for matrix expansion
        - exps: nilpotency class list for EpimorphismPGroup approximation
        - LogObj: an object provided so that this method can save data as members
    """

    if G not in Groups():
        raise ValueError("G must be a group")
    if hasattr(phi, "parent") and phi.parent() == Hom(G, Ints):
        phi = [gap_ZZ_to_int(phi(g)) for g in G.gens()]
    elif not isinstance(phi, list) and not isinstance(phi, tuple):
        raise ValueError("phi must be a morphism G --> Z")
    if abs(gcd(phi)) != 1:
        raise ValueError("phi must be surjective")

    x = lift
    LogObj.lift = x
    Log(INFO, f"Lift: {x}")
    minus_phi = [-z for z in phi]

    B = M.base_ring()
    F = B.group()
    from sage.groups.indexed_free_group import IndexedFreeGroup
    if not isinstance(F, IndexedFreeGroup):
        raise ValueError(f"The base group of M must be an IndexedFreeGroup (got {F})")

    k = M.nrows()

    A0, bias0 = matrix_expansion(G, phi, lift, M, n, LogObj, uniform_bias=True)
    A1, bias1 = matrix_expansion(G, minus_phi, lift**-1, M, n, LogObj, uniform_bias=True)

    if opt.ASK_INPUT and input("Keep a copy of the expanded matrices? "):
        LogObj.expanded_matrix0 = copy(A0)
        LogObj.expanded_matrix1 = copy(A1)

    ord0 = rk * (n + bias0) - von_neumann_rank(G, A0, exps, LogObj)
    ord1 = rk * (n + bias1) - von_neumann_rank(G, A1, exps, LogObj)

    Log(INFO, f"Valuations: {rational_n(ord0)}, {rational_n(ord1)}")
    LogObj.valuation = (ord0, ord1)

    deg = -(ord0 + ord1)
    Log(INFO, f"Degree: {rational_n(deg)}")
    return deg
