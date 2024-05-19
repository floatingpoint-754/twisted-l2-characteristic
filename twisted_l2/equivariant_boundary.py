"""
Provides various utilities for the construction of Z[G]-chain complexes.
"""

from sage.all import *
from random import randrange
from sage.groups.finitely_presented import wrap_FpGroup
import regina
from .utils import to_indexed_free_group, tietze_to_syllables
from .logger import Log, SILENT, INFO, DEBUG

__all__ = [
    "gap_member",
    "ComponentObjectWrapper",
    "face_lattice_to_hap",
    #"regina_rc_to_face_lattice",
    #"verify_face_lattice",
    "regina_tri_to_face_lattice",
    "surface_face_lattice",
    "cw_complex",
    "simplify_cw",
    "equivariant_cc",
    "face_lattice_to_cc",
    "get_fundamental_group",
    "get_differentials",
    #"extract_element",
    #"extract_all_elements",
    "boundary_operator",
    "boundary_operator_row",
    "isosig_RT5",
]

gap_member = libgap.eval("function(obj, str) return obj!.(str); end;")

class ComponentObjectWrapper:
    """
    A class that turns a GAP component object into a Python object with members.
    
    This is needed because GAP component objects, while similar to records,
    cannot be accessed from Sage with the bracket [] operator.
    This is NOT done recursively: any nested component objects remain GAP objects.
    """
    def __init__(self, comp_obj):
        self.__dict__.update({s: gap_member(comp_obj,s) for s in comp_obj.NamesOfComponents()})
        
    def __repr__(self):
        s = "\n ".join(f"{k} : {repr(v)}" for k, v in self.__dict__.items())
        return f"<{s}>"

def face_lattice_to_hap(fl):
    """
    Converts a 0-indexed face lattice to HAP 1-indexed notation.
    """
    dim = len(fl) - 1
    FL = [[[1,0] for _ in fl[0]]]
    for d in range(1, dim+1):
        FL.append([[len(bnd)] + [x+1 for x in bnd] for bnd in fl[d]])
    FL.append([])
    return FL

def regina_rc_to_face_lattice(T, hap = True):
    """
    Returns the face lattice associated to a Regina triangulation.
    
    Arguments:
    
        - T: a Regina Triangulation object. 
        - hap: use the HAP 1-indexed format.
        
    Conditions:
    
        This function assumes that "T" is a regular complex
        (i.e. no simplex has repeated vertices).
        
    Returns:
    
        A list containing, for each dimension,
            a list containing, for each cell C of that dimension,
                a list of integers: the indices of the cells in the boundary of C.
        If "hap" is True, then the indices start from 1, and each list
        of indices is of the form [k, a_1, a_2, ..., a_k],
        where the a_i are cells in the boundary.
    """
    n = T.dimension
    fv = T.fVector()
    if hap:
        FL = [[[1,0] for i in range(fv[0])]]
        for d in range(1,n):
            FL.append([[d+1] + [spx.face(d-1, i).index()+1 for i in range(d+1)] for spx in T.faces(d)])
        FL.append([[n+1] + [T.simplex(k).face(n-1, i).index()+1 for i in range(n+1)] for k in range(fv[n])])
        FL.append([])
        return FL
    else:
        FL = [[[] for i in range(fv[0])]]
        for d in range(1,n):
            FL.append([[spx.face(d-1, i).index() for i in range(d+1)] for spx in T.faces(d)])
        FL.append([[T.simplex(k).face(n-1, i).index() for i in range(n+1)] for k in range(fv[n])])
        return FL

def verify_face_lattice(FL):
    """Verifies that the face lattice returned by regina_rc_to_face_lattice represents a regular CW complex."""
    for fli in FL:
        for lst in fli:
            if len(set(lst[1:])) != lst[0]:
                # print(lst)
                return False
    return True

def regina_tri_to_face_lattice(T, ideal=True, simplify = True):
    """
    Returns the face lattice associated to a Regina triangulation.
    
    Arguments:
    
        - T: a Regina Triangulation object. 
        - ideal: truncate ideal vertices.
        - simplify: simplify the triangulation (uses intelligentSimplify() from Regina).
        
    Conditions:
    
        The dimension of T must be at most 4.
        
    Returns:
    
        A face lattice, that is:
        a list containing, for each dimension,
            a list containing, for each cell C of that dimension,
                a list of integers: the indices of the cells in the boundary of C.
        Each list of indices is of the form [k, a_1, a_2, ..., a_k],
        where the a_i are cells in the boundary (indices start from 1).
        
    Exceptions:
    
        - ValueError: if the dimension of T is greater than 4.
    """
    n = T.dimension
    if n >= 5:
        raise ValueError("Dimension of T must be at most 4")
    if ideal:
        T.idealToFinite()
    if simplify:
        T.intelligentSimplify()
    T.barycentricSubdivision()
    return regina_rc_to_face_lattice(T)

def surface_face_lattice(g, hap=True):
    """Returns the face lattice associated to a closed surface of genus g."""
    if g <= 0:
        raise ValueError("Genus must be positive")
    # 0-cells: 0 (vertex), 1 (center), 2..2g+1 (edges)
    # 1-cells: 0..4g-1 (inner), 4g..8g-1 (outer)
    # 2-cells: 0..4g-1
    
    # 0-indexed format
    FL0 = [[] for i in range(2*g+2)]
    FL1 = [t for i in range(1,g+1) for t in ([0,2*i],[0,2*i+1],[0,2*i],[0,2*i+1])]
    FL1.extend(t for i in range(1,g+1) for t in ([1,2*i],[1,2*i],[1,2*i+1],[1,2*i+1]))
    FL2 = [t for i in range(g) for t in ([4*i+0, 4*i+1, 4*(g+i)+1, 4*(g+i)+2],
                                         [4*i+1, 4*i+2, 4*(g+i)+1, 4*(g+i)+3],
                                         [4*i+2, 4*i+3, 4*(g+i)+0, 4*(g+i)+3],
                                         [4*i+3, 4*i+4, 4*(g+i)+2, 4*(g+i)+4])]
    FL2[-1] = [4*g-1, 0, 8*g-2, 4*g]
    
    if not hap:
        return [FL0, FL1, FL2]
    
    fl0 = [[1,0] for _ in FL0]
    fl1 = [[len(x)]+[i+1 for i in x] for x in FL1]
    fl2 = [[len(x)]+[i+1 for i in x] for x in FL2]
    return [fl0, fl1, fl2, []]
    
def cw_complex(fl):
    return libgap.RegularCWComplex(fl)

def simplify_cw(cw):
    return libgap.SimplifiedRegularCWComplex(cw)
    
def equivariant_cc(cw, gap=False):
    cc = libgap.ChainComplexOfUniversalCover(cw)
    return cc if gap else ComponentObjectWrapper(cc)
    
def face_lattice_to_cc(fl, simplify=True, gap=False):
    cw = cw_complex(fl)
    if simplify:
        cw = simplify_cw(cw)
    return equivariant_cc(cw, gap)
    
def get_fundamental_group(ch):
    return wrap_FpGroup(ch.group)
    
def get_differentials(ch):
    G = get_fundamental_group(ch)
    prop = ch.properties
    dim = 0
    for k,v in prop:
        if k == "dimension": # for complexes
            dim = v.sage()
            break
        if k == "length": # for resolutions
            dim = v.sage() - 1
            break
            
    cc = [boundary_operator(ch, i, G) for i in range(1, dim+1)]
    return cc #to_indexed_free_group(G, cc)
            

# ====================================================================

def extract_element_word(ch, i):
    """Returns a list representing element i (1-based) in the free group."""
    return ch.elts[i-1].UnderlyingElement().LetterRepAssocWord().sage()

def extract_all_element_words(ch):
    """Returns a list of lists representing all group elements in the free group."""
    return [x.UnderlyingElement().LetterRepAssocWord().sage() for x in ch.elts]

def boundary_operator(ch, n, grp = None, ring = QQ, ask = False):
    """
    Returns the boundary operator between dimensions n and n-1
    of the Z[G] (or Q[G]) chain complex associated to 'name'.
    The map is given as a matrix acting by right multiplication.
    """
    G = get_fundamental_group(ch) if grp is None else grp
    ng = G.ngens()
    F = Groups().free(IntegerRange(ng))
    A = F.algebra(ring)
    h = ch.dimension(n).sage()    # source
    w = ch.dimension(n-1).sage()  # target
    M = matrix(A, h, w)
    bnd = []
    for i in range(h):
        if ask:
            if input(f"Next row? (Y/n) ") == "n":
                bnd.append([])
                continue
        bnd.append(ch.boundary(n, i+1).sage())
        Log(INFO, f"Extracted boundary row #{i} of length {len(bnd[-1])}")
                
    Log(INFO, "Extracting all group elements...")
    elts = extract_all_element_words(ch)
    Log(INFO, f"Extracted all group elements (amount = {len(elts)})")
    for i in range(h):
        for l in bnd[i]:
            cell = l[0]
            elt = l[1]
            a = A.monomial(F(tietze_to_syllables(elts[elt-1])))
            if (cell < 0):
                a = -a
            M.add_to_entry(i, abs(cell)-1, a)
    return M

def boundary_operator_row(ch, n, i, grp = None, ring = QQ):
    """
    Returns one specific row of the boundary operator.
    For debug purposes, when computation is so intensive that
    we want to know where it freezes.
    """
    G = get_fundamental_group(ch) if grp is None else grp
    ng = G.ngens()
    F = Groups().free(IntegerRange(ng))
    A = F.algebra(ring)
    h = ch.dimension(n).sage()    # source
    w = ch.dimension(n-1).sage()  # target
    M = matrix(A, 1, w)
    
    bnd = ch.boundary(n, i+1)
    Log(DEBUG, "GAP call completed")
    bnd = bnd.sage()
    Log(INFO, f"Extracted boundary #{i} of length {len(bnd[-1])}")
                
    elts = {l[1]-1: ch.elts[l[1]-1].UnderlyingElement().LetterRepAssocWord().sage() for l in bnd}
        #extract_element(ch, l[1]) for l in bnd}
    Log(INFO, f"Extracted all group elements (amount = {len(elts)})")
    for l in bnd:
        cell = l[0]
        elt = l[1]
        a = A.monomial(F(tietze_to_syllables(elts[elt-1])))
        if (cell < 0):
            a = -a
        #print(i, abs(cell)-1, a)
        M.add_to_entry(0, abs(cell)-1, a)
    return M

# isomorphism signature of ideal triangulation of RT5 manifold
isosig_RT5 = (
    "-cadvvAMzzLzvAvvLvvvzvzLPQvzLvzMLvQvALvLQQMQAwPwwLAQwMQQwMPQAQMw"
    "QPwQMAQPMAQQAMQvwPMzMvQPMPvPQPvvQPQAzvAMwQQvQQMQAQQQLwQQvzQwwLQw"
    "QQQLQzQQQAzQQAQAQQQQQQPQLMQQQQQPQQQPQQQQQQLQQQQQMQQPQQQQQAQQPQPQ"
    "QQQQcadafadafafakajaoajazaIaOavaFauasataGaEaZa3atauauaTaUaVaEa9a"
    "La0a5acbfbjb5aBaCaDadbCaOaDaLaDadbeb9aMa1aIaJaKahbJaKaKa8albhbgb"
    "ibsbtbubBb2aCbzbFbQaRaqbSaRa3aSa2aSawbabUaVaAbVaxbXaYaZaLbYaZa+a"
    "ZaybsblbmbkbvbxbwbDbBbEbQbRbPb5a6a7a6a7aqb7atbdbebZbwbHbzbHbdbnb"
    "ObMbnbWbJbUbebCbRbhbAbib+biblbmbLbmbZbdcgc9b3bKbhcGbncWbVbPbVbTb"
    "tb+bububwbxbZbYbxbvcEb-b+b2b3b4bwcvcEbvcDbacRbbcEbdclcecpcocVbqc"
    "ObtcCc9b8bHcUbkcNc7bWbObDcGcNc7bWbscTbtczcFcbcacQcicocsctcTc9b8b"
    "Vb7b7bcc-b0bAcyc5bucVccc-bUcDcucLc6btc3b4bUc4b5bTcScWcycAcXcLcuc"
    "9brcTc9bQcAcccUcTcbcocUcvc1c0cBc0cCclcgcHchcicGc0chcicicscrcScmc"
    "EcNcJcKcIcRcBcCcEcScrcOcRc4cscOcpcqcHcqcCcNcLcEcSctcLcXc2cxcMcMc"
    "2c6cZcWcAc5cFcRcQc0cWcCcYcOcPcGcScKcQcRcOcPcKc8cNcJcKcMcYcKcWcXc"
    "6c2cPcYc6cRc3c7c1cZc5c9c9c+c9c7c1c5c7c5c6c-c-c8c+c4c-c8c+c7c9c9c"
    "+c-c-c4bieaiyaqboagfgfgfgfaagfaaaagfIe0iOigfgfaaaa4baiqb2agaca4b"
    "Yaieai4haaaaaa4h4bieaiqiyaueqbueoaUc+bOfyaqb4bieai2fyaqboaec4dGe"
    "oekeEgEaAaGhoaaagfaa4bieWgaiyaogqbogoaQfIeyaqbgfoaSh4bieggaaya2a"
    "OigagfSh2igi+h8gwbobmgCimaaaaaaa4bieaiyaqbOeoagfqiobaagfQiaaQi6e"
    "6b6b+h6baaggggieue4h4bafaiaaqb4biecayaYaaaiaaaOisiaaieaaWfWfWfWf"
    "ai4b+baiqb4bieOfOfyaaaYg2f2f4byaqbweweqeweGaiaogkeiekeAg6hqcaaei"
    "6hIhSeaaoeEgaaoakbaaGgogIhaamgaacagdobgdad6hmdgioeaaycmaidweaa+h"
    "8g4eqeqbCiaa4bUf+baigfaa4d+b6bQiShSeqb6bieai8doaefadad2iOhOh8gke"
    "Ci6b4gib2aibaiKgoaibieYaOfIhiaiaEciaUeYa4bggggaiggaa2aqbcaubub2g"
    "ygogca+aGaKaoaibafafIhoaygoaaaubyg4baiqeqbqeygygIeMhyaIh+aqcaiWc"
    "0iQi+h2f2dyawe0iafUeag2agaef4bieoa+dGa+dWf0cueGaaa+dyaqb+d2doawi"
    "+gobkeya+gWgcaycmdaaob+baa6boe8dgaGaod2gecieiaqcMe0cYaGhqecdsiMh"
    "Wi2gUf2c"
)










































































































































