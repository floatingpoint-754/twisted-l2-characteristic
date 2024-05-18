"""
Defines functions computing the twisted rank of matrices with coefficients in a group ring.
"""

from sage.all import *
from sage.groups.group import is_Group
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
    
    if not is_Group(G):
        raise ValueError("G must be a group")
    
    
    B = M.base_ring()
    Fb = B.group()
    from sage.groups.indexed_free_group import IndexedFreeGroup
    
    k = M.nrows()
    h = M.ncols()
    
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
    
def determinant_degree(G, phi, lift, M, n, exps, LogObj = lambda: None):
    """
    Computes the degree of the Dieudonné determinant of a square matrix with coefficients in the group ring Q[G].
    
    Arguments:
    
        - G: finitely presented group
        - phi: surjective morphism G --> Z (can be given as a list of integers of length G.ngens())
        - M: square matrix over Q[G] (given as a matrix over the free group)
        - n: upper bound for matrix expansion
        - exps: nilpotency class list for EpimorphismPGroup approximation
    
    Complexity:
    
        The matrix M is enlarged by a factor linear in n and possibly exponential in p and c.
    """
    
    if not is_Group(G):
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
    
    def val_F(f):
        return sum((s*phi[j] for j,s in f.to_word_list()), 0)
    
    # Now we need to find dim_D(K) coker(M), that is deg_K Det_D(G)(M)
    # We use Taihei Oki's matrix expansion algorithm (arXiv:1907.04512, 2019)
    # First, we multiply the rows on the left by adequate multiples of x,
    # such that all of M's entries have positive valuation
    
    valmin = lambda y: min((val_F(f) for f,q in y), default=0)
    minrows = [min((valmin(M[i,j]) for j in range(k)), default=0) for i in range(k)]
    Log(DEBUG, f"Minimum exponent of each row: {minrows}, total: {sum(minrows)}")
    
    for i in range(k):
        for j in range(k):
            M[i,j] = B.monomial(x**(-minrows[i]))*M[i,j]
                
    valmax = lambda y : max((val_F(f) for f,q in y), default=0)
    
    maxval = max((max((valmax(M[i,j]) for j in range(k)), default=0) for i in range(k)), default=0)
    mu = min(n, maxval*k)
    
    Log(INFO, f"Maximum valuation of entries = {maxval}, should expand matrix by {maxval*k}, expanding by {mu} instead")
    if mu == 0:
        Log(INFO, "Valuation: 0 (skipping this...)")
        return sum(minrows)*-2
    
    Ax = matrix(B, mu*k, mu*k)
    
    ProgressBar(DEBUG, PB_START)
    for i in range(k):
        for j in range(k):
            l = [B.monomial(x**s) * M[i,j] for s in range(mu)]
            for s in range(mu):
                for h, q in l[s]:
                    H = val_F(h)
                    if H in range(mu):
                        Ax.add_to_entry(i+k*s, j+k*H, q * B.monomial(h * x**(-H)))
            
            ProgressBar(DEBUG, (i*k+j)/(k*k))
    
    ProgressBar(DEBUG, PB_STOP)
    Log(DEBUG, "expanded matrix done")
    
    if opt.ASK_INPUT and input("Keep a copy of the expanded matrix? "):
        LogObj.expanded_matrix = copy(Ax)
        
        
    rnk = von_neumann_rank(G, Ax, exps, LogObj)
    valn = mu*k - rnk
    Log(INFO, f"Valuation: {rational_n(valn)}")
    
    LogObj.valuation = valn
    
    deg = (valn + sum(minrows))*-2
    Log(INFO, f"Degree: {rational_n(deg)}")
    return deg
