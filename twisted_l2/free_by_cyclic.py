"""
Provides the FreeByCyclic class, used to make classifying spaces for free-by-cyclic groups.
"""

from sage.all import *

__all__ = ["FreeByCyclic"]

GEN_EDGE_BACK = 0
GEN_EDGE_FRONT = 1
GEN_EDGE_TO_C = 2
GEN_EDGE_CBACK = 3
GEN_EDGE_CFRONT = 4
GEN_FACE_QBACK = 0
GEN_FACE_QFRONT = 1
GEN_FACE_TBACK = 2
GEN_FACE_TFRONT = 3

class FreeByCyclic:
    """Creates a face lattice which is a classifying space for a free-by-cyclic group."""
    def __init__(self, F, phi):
        """
        Constructor for FreeByCyclic class.
        
        Arguments:
        
            - F: free group.
            - phi: endomorphism of F.
        """
        self.F = F
        self.N = F.ngens()
        self.phi = phi
        self.rels = [phi(g).Tietze() for g in F.gens()]
        self.clen = [0]
        for r in self.rels:
            self.clen.append(self.clen[-1] + len(r))
        self.FL = []
    
    def vtx_base(self):
        return 0
    def vtx_midT(self):
        return 1
    def vtx_midG(self,n):
        return 2*n+2
    def vtx_faceG(self,n):
        return 2*n+3
    
    def edge_backT(self):
        return 0
    def edge_frontT(self):
        return 1
    def edge_G(self, n, edgetype):
        return 5*n+2+edgetype
    def edge_R(self, n, k):
        return 5*self.N+2 + 2*self.clen[n]+n + k
    def edge_RL(self, n):
        return self.edge_R(n+1,0)-1
    
    def face_G(self, n, facetype):
        return 4*n+facetype
    def face_R(self, n, k):
        return 4*self.N + 2*self.clen[n] + k
    
    def make_verts(self):
        self.FL.append([[] for _ in range(2*self.N + 2)])
        
    def make_edges(self):
        l = [[] for _ in range(5*self.N+2 + 2*self.clen[self.N]+self.N)]
        l[self.edge_backT()] = [self.vtx_base(), self.vtx_midT()]
        l[self.edge_frontT()] = [self.vtx_base(), self.vtx_midT()]
        
        for i in range(self.N):
            l[self.edge_G(i, GEN_EDGE_BACK)] = [self.vtx_base(), self.vtx_midG(i)]
            l[self.edge_G(i, GEN_EDGE_FRONT)] = [self.vtx_base(), self.vtx_midG(i)]
            l[self.edge_G(i, GEN_EDGE_TO_C)] = [self.vtx_faceG(i), self.vtx_midG(i)]
            l[self.edge_G(i, GEN_EDGE_CBACK)] = [self.vtx_faceG(i), self.vtx_midT()]
            l[self.edge_G(i, GEN_EDGE_CFRONT)] = [self.vtx_faceG(i), self.vtx_midT()]
            
        for i, r in enumerate(self.rels):
            for k in range(len(r)+1):
                l[self.edge_R(i, 2*k)] = [self.vtx_faceG(i), self.vtx_base()]
            for k, g in enumerate(r):
                l[self.edge_R(i, 2*k+1)] = [self.vtx_faceG(i), self.vtx_midG(abs(g)-1)]
                
        self.FL.append(l)
        
        
    def make_faces(self):
        l = [[] for _ in range(4*self.N + 2*self.clen[self.N])]
        
        for i in range(self.N):
            l[self.face_G(i, GEN_FACE_QBACK)] = [self.edge_G(i, GEN_EDGE_BACK),
                                                self.edge_G(i, GEN_EDGE_CBACK),
                                                self.edge_G(i, GEN_EDGE_TO_C),
                                                self.edge_backT()]
            l[self.face_G(i, GEN_FACE_QFRONT)]= [self.edge_G(i, GEN_EDGE_FRONT),
                                                self.edge_G(i, GEN_EDGE_CFRONT),
                                                self.edge_G(i, GEN_EDGE_TO_C),
                                                self.edge_backT()]
            l[self.face_G(i, GEN_FACE_TBACK)]= [self.edge_G(i, GEN_EDGE_CBACK),
                                                self.edge_R(i,0),
                                                self.edge_frontT()]
            l[self.face_G(i, GEN_FACE_TFRONT)]= [self.edge_G(i, GEN_EDGE_CFRONT),
                                                self.edge_RL(i),
                                                self.edge_frontT()]
            
        for i, r in enumerate(self.rels):
            for k, g in enumerate(r):
                edges = [self.edge_G(abs(g)-1, GEN_EDGE_BACK), self.edge_G(abs(g)-1, GEN_EDGE_FRONT)]
                if g < 0:
                    edges[0], edges[1] = edges[1], edges[0]
                    
                l[self.face_R(i, 2*k)] = [self.edge_R(i, 2*k), self.edge_R(i, 2*k+1), edges[0]]
                l[self.face_R(i, 2*k+1)] = [self.edge_R(i, 2*k+1), self.edge_R(i, 2*k+2), edges[1]]
                
        self.FL.append(l)
    
    def make_fl(self):
        """
        Returns a face lattice for the classifying space of this free-by-cyclic group.
        """
        self.make_verts()
        self.make_edges()
        self.make_faces()
        return self.FL
        
    def _test_me(self):
        """
        This function is a work in progress.
        
        Code needs to be refactored/updated; this doubles as a test for the RegularCWComplex class.
        """
        fl = self.make_fl()
        
        import regular_cw_complex as rcc
        import equivariant_boundary as eqb
        
        rc = rcc.RegularCWComplex(fl)
        FL, _ = rc.face_lattice_orientation()
        
        print("Computed face lattice:", FL)
        
        rc = eqb.cw_complex(FL)
        print(f"Size(rc) = {libgap.Size(rc)}, Size(rcc) = {libgap.Size(libgap.ContractedComplex(rc))}")
        
        ch_gap = eqb.equivariant_cc(rc)
        ch = eqb.ComponentObjectWrapper(ch_gap)
        G = eqb.get_fundamental_group(ch)
        
        print(G)
        print("Abelian invariants:", G.abelian_invariants())
        print("Dimensions:", [ch.dimension(j) for j in range(3)])
        Ch = [eqb.boundary_operator(ch, 1, grp=G), eqb.boundary_operator(ch, 2, grp=G)]
        
        # make some loops
        T = [self.edge_backT()+1, self.edge_frontT()+1]
        Gs = [[self.edge_G(i, GEN_EDGE_BACK)+1, self.edge_G(i, GEN_EDGE_FRONT)+1] for i in range(self.N)]
        
        edges = [T]+Gs
        edge_to_word = eqb.ComponentObjectWrapper(ch.group).edgeToWord
        loops = [G(libgap.LetterRepAssocWord(libgap.UnderlyingElement(edge_to_word(e[0])*edge_to_word(e[1])))) for e in edges]
        print("Loops =", loops)
        
        return {"group": G, "complex": Ch, "loops": loops}
        
    @staticmethod
    def random_automorphism(n, length=4, seed=None, print_latex=True):
        """
        Generates a random automorphism of F_n.
        
        The automorphism is returned as a product of generators
        of the form tau_i, sigma_ij, eta_ij.
        
        Arguments:
        
            - n: rank of the free group.
            - length: the average length of the automorphism.
            - seed: an integer, can be specified to make the output deterministic.
            - print_latex: print the automorphism in LaTeX notation.
        """
        import random
        if seed is None:
            import sys
            seed = random.randrange(sys.maxsize)
        rng = random.Random()
        rng.seed(seed, version=2)
        FreeByCyclic.random_automorphism.seed = seed
        
        F = FreeGroup(n)
        gens = list(F.gens())
        
        def tau(i):
            l = gens[:]
            l[i] = l[i]**-1
            return (f"\\tau_{{{i+1}}}"), F.hom(l)
        
        def sigma(i,j):
            l = gens[:]
            l[i],l[j] = l[j],l[i]
            return (f"\\sigma_{{{i+1},{j+1}}}"), F.hom(l)
        
        def eta(i,j):
            l = gens[:]
            l[j] = l[j]**-1
            l[i] = l[j]*l[i]
            return (f"\\eta_{{{i+1},{j+1}}}"), F.hom(l)
            
        gl = [tau(i) for i in range(n)]
        gl.extend(sigma(i,j) for i in range(n) for j in range(n) if i < j)
        gl.extend(eta(i,j) for i in range(n) for j in range(n) if i != j)
        
        none_count = n * (3*n-1)/(2*length)
        if none_count < 1:
            rep = int(1/none_count + 1/2)
            gl = gl*rep
            none_count = 1
        else:
            none_count = int(none_count + 1/2)
            
        gl.extend(None for _ in range(none_count))
        
        aut = F.hom(gens)
        inv = F.hom(gens)
        while True:
            g = rng.choice(gl)
            if g is None:
                break
            aut = aut*g[1]
            inv = g[1]*inv
            if print_latex:
                print(g[0], end=" ")
        if print_latex:
            print()
        
        return F, F.hom([aut(x) for x in gens]), F.hom([inv(x) for x in gens])
