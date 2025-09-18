"""
Provides methods for the computation of the twisted L^2-Euler characteristic.
"""

from sage.all import *
from . import twisted_rank as trk
from .utils import tietze_to_syllables, ga_adjoint, Ints, gap_ZZ_to_int, rational_n
from .logger import Log, SILENT, MINIMAL, INFO, DEBUG
from . import configs

__all__ = [
    "make_phi_from_coordinates",
    "laplacians",
    "twisted_l2_characteristic",
    "characteristic",
    "twisted_l2_betti_number",
    "betti_number",
    "twisted_l2_betti_number_kth",
    "betti_number_kth",
    "get_twisted_l2_logs",
    "draw_seminorm"
]

class LogObject:
    """An object to which members can be added at will."""
    pass

def make_phi_from_coordinates(G, v, *, as_list=False):
    """
    Makes a group homomorphism phi: G --> Z given by the list of integers v.
    
    Arguments:
    
        - G: finitely presented group.
        - v: a list or tuple of integers, the coordinates describing phi.
        
    Keyword arguments:
    
        - as_list:
        
            - if True, return the images of the generators of G under phi.
            - if False, return phi as a Sage morphism.
    
    Exceptions:
    
        - ValueError: if len(v) is not equal to rk(G).        
    """
    
    r = G.abelian_invariants().count(0)
    if len(v) != r:
        raise ValueError("The length of v must be equal to the rank of G")
        
    abhom = libgap.MaximalAbelianQuotient(G)
    ab = abhom.Image()    # GAP gives a FP group where torsion elements come first
    n = ab.GeneratorsOfGroup().Size().sage()
    abphi = libgap.GroupHomomorphismByImages(ab, Ints, ab.GeneratorsOfGroup(), [Ints([0])]*(n-r)+[Ints([i]) for i in v])
    phigap = libgap.CompositionMapping(abphi, abhom)
    gens_images = [phigap.Image(g) for g in G.gens()]
    
    if as_list:
        return [gap_ZZ_to_int(Ints(m)) for m in gens_images]
    else:
        return G.Hom(Ints)(gens_images)
    
def laplacians(cc):
    """
    Returns a list of the combinatorial Laplacian operators of the chain complex cc.
    
    Arguments:
    
        - cc: a chain complex (list) of matrices over Q[F], where F is a free group.
    
    Returns:
    
        A list of length len(cc)+1 containing square matrices, as follows:
            [cc[0]' cc[0], cc[1]' cc[1] + cc[0] cc[0]', ..., cc[n] cc[n]']
    """
    
    d = len(cc)-1
    laps = [None for _ in range(d+2)]
    Ct = [ga_adjoint(M) for M in cc]
    laps[0] = Ct[0] * cc[0]
    laps[-1] = cc[-1] * Ct[-1]
    for i in range(1, d+1):
        laps[i] = Ct[i]*cc[i] + cc[i-1]*Ct[i-1]
    return laps
        
    
def twisted_l2_characteristic(G, v, cc, n, exps):
    """
    Computes the twisted L^2-Euler characteristic of a chain complex with coefficients in the group ring Q[G].
    
    Arguments:
    
        - G: finitely presented group.
        - v: coordinates of a morphism phi: G --> Z in H^1(G; Z).
        - cc: chain complex over Q[G] given as a list of matrices
            (cc[i] is the matrix between dimensions (i+1) and (i)).
        - n: multiplier for approximation by a matrix in Q[K]
            (K = ker phi).
        - exps: nilpotency class list for EpimorphismPGroup approximation.
        
    Complexity:
    
        Each matrix is enlarged by a factor linear in n and possibly exponential in p and c.
    """
    if len(cc) == 0:
        Log(MINIMAL, f"Characteristic        : 0")
        return 0
    
    fac = gcd(v)
    if fac == 0:
        Log(MINIMAL, f"Characteristic        : 0")
        return 0
    
    phi = make_phi_from_coordinates(G, [vi/fac for vi in v]) # make a primitive phi
    x = phi.lift(Ints([1]))
    phi = [gap_ZZ_to_int(phi(g)) for g in G.gens()]
    d = len(cc)
    for i in range(d-1):
        if cc[i].nrows() != cc[i+1].ncols():
            raise ValueError("cc has inconsistent dimensions")
    
    if not isinstance(exps, tuple):
        exps = tuple(exps)
        
    # convert cc to IndexedFreeGroup
    F = Groups().free(IntegerRange(G.ngens()))
    QF = F.algebra(QQ)
    if cc[0].base_ring() is not QF:
        Log(DEBUG, "Converting to IndexedFreeGroup...")
        cc = to_indexed_free_group(G, cc)
        Log(DEBUG, "Done")
        
    lift = F(tietze_to_syllables(x.Tietze()))
    Log(DEBUG, "Computing Laplacians...")
    
    Laps = laplacians(cc)
    Log(DEBUG, f"Dimensions of the Laplacians: {[M.dimensions() for M in Laps]}")
        
    logs = [LogObject() for L in Laps]
    laplacian_degs = []
    
    for j, L in enumerate(Laps):
        Log(INFO, f" Dimension {j} ".center(35, "="))
        Log(INFO, "")
        det_deg = trk.determinant_degree(G, phi, lift, L, n, exps, LogObj = logs[j])
        laplacian_degs.append(det_deg)
        Log(INFO, "")
    
    for l, ld in zip(logs, laplacian_degs):
        l.degree = ld
        
    #Log(MINIMAL, f"GCD of entries of v: {fac}")
    Log(MINIMAL, f"Degrees of Laplacians : {[rational_n(ld) for ld in laplacian_degs]}")
    Log(MINIMAL, f"Valuations            : {[rational_n(l.valuation) for l in logs]}")
    Log(MINIMAL, f"Quotient sizes        : {[l.Lsize for l in logs]}")
    
    terms = [i*x for i,x in enumerate(laplacian_degs)]
    twisted_l2_characteristic.logs = logs
    result = -Integer(fac)/2 * (sum(terms[::2]) - sum(terms[1::2]))
    Log(MINIMAL, f"Characteristic        : {rational_n(result)}")
    return result

characteristic = twisted_l2_characteristic

def get_twisted_l2_logs():
    """
    Returns a list of "log" objects associated to the last run of twisted_l2_characteristic,
    or None if the function has not run yet.
    
    The members of these objects contain information about each run of determinant_degree.
    """
    if hasattr(twisted_l2_characteristic, "logs"):
        return twisted_l2_characteristic.logs
    else:
        return None

def twisted_l2_betti_number(G, v, cc, k, n, exps):
    """
    Computes the i-th twisted L^2-Betti number of a chain complex with coefficients in the group ring Q[G].

    Arguments:

        - G: finitely presented group.
        - v: coordinates of a morphism phi: G --> Z in H^1(G; Z).
        - cc: L^2-acyclic chain complex over Q[G] given as a list of matrices
            (cc[i] is the matrix between dimensions (i+1) and (i)).
        - k: dimension in which to calculate the Betti number.
        - n: multiplier for approximation by a matrix in Q[K]
            (K = ker phi).
        - exps: nilpotency class list for EpimorphismPGroup approximation.

    Complexity:

        Polynomial in the size of cc[i], polynomial in n, and possibly exponential in p and c.
    """
    if k >= len(cc) or k < 0:
        Log(MINIMAL, f"b_{k}          : 0")
        return 0

    fac = gcd(v)
    if fac == 0:
        Log(MINIMAL, f"b_{k}          : 0")
        return 0

    phi = make_phi_from_coordinates(G, [vi/fac for vi in v]) # make a primitive phi
    x = phi.lift(Ints([1]))
    phi = [gap_ZZ_to_int(phi(g)) for g in G.gens()]
    d = len(cc)

    for i in range(d-1):
        if cc[i].nrows() != cc[i+1].ncols():
            raise ValueError("cc has inconsistent dimensions")

    if not isinstance(exps, tuple):
        exps = tuple(exps)

    # convert cc to IndexedFreeGroup
    F = Groups().free(IntegerRange(G.ngens()))
    QF = F.algebra(QQ)
    if cc[0].base_ring() is not QF:
        Log(DEBUG, "Converting to IndexedFreeGroup...")
        cc = to_indexed_free_group(G, cc)
        Log(DEBUG, "Done")

    lift = F(tietze_to_syllables(x.Tietze()))

    rk = sum(cc[j].ncols() * (-1)**(k-j) for j in range(k+1))
    Log(MINIMAL, f"Rank : {rk}")

    logobj = LogObject()
    minors = trk.determinant_degree_minors(G, phi, lift, cc[k], rk, n, exps, LogObj = logobj)

    Log(MINIMAL, f"Best minors    : {logobj.best_minors[0]} (rows), {logobj.best_minors[1]} (columns)")
    Log(MINIMAL, f"Valuation      : {rational_n(logobj.valuation)}")
    Log(MINIMAL, f"Quotient sizes : {logobj.Lsize}")
    Log(INFO,    f"Degrees        : {minors}")

    result = minors[-1]
    Log(MINIMAL, f"b_{k}          : {rational_n(result)}")
    return result

betti_number = twisted_l2_betti_number

def twisted_l2_betti_number_kth(k):
    """
    Returns a curried version of `twisted_l2_betti_number` with the same signature as `twisted_l2_characteristic`,
    computing the k-th twisted L^2-Betti number.
    """
    return lambda G, v, cc, n, exps: twisted_l2_betti_number(G, v, cc, k, n, exps)

betti_number_kth = twisted_l2_betti_number_kth

def draw_seminorm(G, cc, n, exps, rad = 1, seminorm = twisted_l2_characteristic):
    """
    Constructs a matrix of values of the twisted L^2 Euler characteristic (or another seminorm) at integer points.
    
    Arguments:
    
        - G: finitely presented group.
        - cc: chain complex over Q[G] given as a list of matrices (cc[i] is the matrix between dimensions (i+1) and (i)).
        - n: multiplier for approximation by a matrix in Q[K] (K = ker phi).
        - exps: nilpotency class list for EpimorphismPGroup approximation.
        - rad: radius of the matrix (e.g. rad=3 returns a 7x7 matrix).
        - seminorm: the seminorm to draw (default: `twisted_l2_characteristic`, may also be `twisted_l2_betti_number_kth(k)`).
    
    Conditions:
    
        Only dim(H^1) <= 2 is supported.
    
    Exceptions:
    
        - NotImplementedError: if rk(G) > 2.
    """
    r = G.abelian_invariants().count(0)
    if r > 2:
        raise NotImplementedError("Rank of G must be 2 or less")
    if r == 0:
        return matrix([[0]])
    
    nr = 1 if r == 1 else 2*rad+1
    nc = 2*rad+1
    M = matrix(RDF, nr, nc)
    for i in range(nr):
        x = i - nr // 2
        for j in range(nc):
            y = -(j - nc // 2)
            M[i,j] = twisted_l2_characteristic(G, [y,x], cc, n, exps)
    
    return M
