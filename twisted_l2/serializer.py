from sage.all import *
from .utils import tietze_to_syllables
import json

__all__ = ["save_to_file", "load_from_file"]

class FpGroupSerializer:
    def __init__(self, G):
        self.ngens = G.ngens()
        self.relations = tuple(r.Tietze() for r in G.relations())
    
    def __str__(self):
        l = [f"Finitely presented group on {self.ngens} generators with relations:"]
        l.extend(str(r) for r in self.relations)
        return "\n\t".join(l)
    
    def group(self):
        F = FreeGroup(self.ngens)
        return F / self.relations

class FreeGroupAlgebraMatrixSerializer:
    def __init__(self, M):
        self.R = M.base_ring().base_ring()
        self.Fn = M.base_ring().group().ngens()
        self.dim = M.dimensions()
        self.entries = [[{g.Tietze():v for g, v in M[i,j]} for j in range(M.ncols())] for i in range(M.nrows())]
    
    def __str__(self):
        return f"{self.dim[0]} x {self.dim[1]} matrix on Algebra of Free Group on {self.Fn} generators over {self.R}"
    
    def matrix(self):
        F = FreeGroup(self.Fn)
        A = F.algebra(self.R)
        return matrix(A, self.dim[0], self.dim[1], lambda i,j: A.sum_of_terms({F(t):v for t, v in self.entries[i][j].items()}, distinct=True))


def fp_group_to_dict(G):
    return {"ngens": int(G.ngens()), "relations": [intify(r.Tietze()) for r in G.relations()]}
    
def dict_to_fp_group(d):    
    F = FreeGroup(d["ngens"])
    return F / d["relations"]

def ga_matrix_to_dict(M):
    A = M.base_ring()
    R = A.base_ring()
    if R == ZZ:
        Rs = "ZZ"
        encode = lambda v: int(v)
    elif R == QQ:
        Rs = "QQ"
        encode = lambda q: (int(q.numerator()), int(q.denominator()))
    else:
        raise ValueError("Group algebra base ring must be ZZ or QQ")
    
    ngens = len(A.group().gens())
    nrows = int(M.nrows())
    ncols = int(M.ncols())
    entries = [[[([int((j+1)*s) for j,s in g.to_word_list()], encode(v)) for g, v in M[i,j]]
                for j in range(M.ncols())]
                for i in range(M.nrows())]
    return {"base_ring": Rs, "ngens": ngens, "nrows": nrows, "ncols": ncols, "entries": entries}
    
def dict_to_ga_matrix(d):
    Rs = d["base_ring"]
    if Rs == "ZZ":
        R = ZZ
        decode = lambda v: Integer(v)
    elif Rs == "QQ":
        R = QQ
        decode = lambda t: Integer(t[0]) / Integer(t[1])
    else:
        raise ValueError("'base_ring' must be 'ZZ' or 'QQ'")
        
    F = Groups().free(IntegerRange(d["ngens"]))
    A = F.algebra(R)
    entries = d["entries"]
    return matrix(A, d["nrows"], d["ncols"], lambda i,j: A.sum_of_terms({F(tietze_to_syllables(t)): decode(v) for t, v in entries[i][j]}, distinct=True))
    
def intify(L):
    if L is None:
        return None
    if not isinstance(L, (list, tuple)):
        return int(L)
    return [intify(l) for l in L]
    
def save_to_file(path, *, group, cc, face_lattice=None):
    gd = fp_group_to_dict(group)
    ch = [ga_matrix_to_dict(m) for m in cc]
    fl = intify(face_lattice)
    d = {"group": gd, "complex": ch}
    if fl is not None:
        d["face_lattice"] = fl
    with open(path, "w") as fil:
        json.dump(d, fil)

def load_from_file(path, as_tuple=True):
    with open(path, "r") as fil:
        d = json.load(fil)
    G = dict_to_fp_group(d["group"])
    cc = [dict_to_ga_matrix(m) for m in d["complex"]]
    
    if as_tuple:
        if "face_lattice" in d:
            return G, cc, d["face_lattice"]
        else:
            return G, cc
    
    D = {"group": G, "complex": cc}
    if "face_lattice" in d:
        D["face_lattice"] = d["face_lattice"]
    return D
        
    
    
    
    
    
    
    
    
    
    
    
    
    
