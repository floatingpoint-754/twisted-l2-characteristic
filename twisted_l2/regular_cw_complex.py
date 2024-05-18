"""
Contains a partial rewriting and extension of HAP's RegularCWComplex type.

This is a WORK IN PROGRESS which may be unneeded in some cases (see 'edgeToWord').
It was originally meant to keep track of particular loops through simplifications.

This file is partially based on HAP source code. The GAP package
HAP is released under GPL 2.0 and maintained by Graham Ellis.
"""

from collections import Counter
from sage.all import *
import copy as cp
    
__all__ = ["RegularCWComplex"]    

class RegularCWComplex:
    """
    A partial rewriting and extension of HAP's RegularCWComplex type.

    This is a WORK IN PROGRESS which may be unneeded in some cases (see 'edgeToWord').
    It was originally meant to keep track of particular loops through simplifications.
    """
    
    def __init__(self, boundaries, orientation = None, coboundaries = None, marked_loops = []):
        """
        Initializes an instance of RegularCWComplex.
        
        Arguments:
        
            - boundaries: a list with one entry for each dimension.
                - boundaries[d]: a list of cells in dimension d.
                    - boundaries[d][i]: a list of (d-1)-cells in the boundary of cell i.
            - orientation: a nested list, where orientation[d][i][k] is the sign
                of cell boundaries[d][i][k] in the boundary of cell i.
                
            OR
            
            - boundaries[d][i]: a dict[int,int]: cell -> orientation.
            - orientation: the string "dict".
        
                
        - marked_loops: a list of loops in the complex.
            - marked_loops[i]: a list containing a 0-cell (or None), then
                1-cells in their order.
        """
        if orientation == "dict":
            self.boundaries = boundaries
            self.oriented = True
        elif orientation is not None:
            self.boundaries = [[{l: orientation[d][i][k]
                for k,l in enumerate(y)}
                for i,y in enumerate(x)]
                for d,x in enumerate(boundaries)]
            self.oriented = True
        else:
            self.boundaries = [[{l: None
                for l in y}
                for y in x]
                for x in boundaries]
            self.oriented = False
        
        if coboundaries is None:
            self.coboundaries = self.calculate_coboundaries()
        else:
            self.coboundaries = coboundaries        
        
        unspooled = False
        if marked_loops:
            if marked_loops[0]:
                if isinstance(marked_loops[0][0], tuple):
                    unspooled = True
            else:
                raise ValueError("Invalid marked_loops")
                
        self.marked_loops = (marked_loops if unspooled
            else [self.unspool_loop(l) for l in marked_loops])
            
        self._NULL_DICT = dict() # used to catch removed cells
    
    def calculate_coboundaries(self):
        cbd = [[set() for _ in x] for x in self.boundaries]
        #print([len(x) for x in self.boundaries])
        #print([len(x) for x in cbd])
        for d,x in enumerate(self.boundaries):
            if d > 0:
                for i,y in enumerate(x):
                    for c in y:
                        #cell c is in boundary of d-cell i
                        #try:
                        cbd[d-1][c].add(i)
                        #except IndexError:
                        #    print("d = %s, len(cbd) = %s, len(cbd[d-1]) = %s, c = %s" % (d, len(cbd), len(cbd[d-1]), c))
                        #    raise IndexError()
        return cbd
    
    def unspool_loop(self, loop):
        """
        Normalizes a loop.
        
        ARGUMENTS:
        - loop: a list in the format [b, e_1, ..., e_k] where b is
            a common 0-cell of e_1 and e_k. Can be None if k > 2.
            
        RETURNS:
        a list of 3-tuples (e,s,t) with
            e: 1-cell edge
            s: start of edge
            t: end of edge
        """
        
        edge_list = []
        start = loop[0]
        if start is None:
            if len(loop) <= 3:
                raise ValueError("loop[0] cannot be None if the path has <= 2 cells")
            for k in self.boundaries[1][loop[1]].keys() & self.boundaries[1][loop[-1]].keys():
                break
            start = k
        for e in loop[1:]:
            for t in self.boundaries[1][e].keys():
                if t != start:
                    break
            edge_list.append((e, start, t))
            start = t
        return edge_list
            
    def ncells(self, d):
        return len(self.boundaries[d])
    
    def dimension(self):
        return len(self.boundaries)-1
        
    def size(self):
        return sum(self.ncells(i) for i in range(self.dimension() + 1))

    def copy(self):
        Y = RegularCWComplex(cp.deepcopy(self.boundaries),
            orientation = "dict",
            coboundaries = cp.deepcopy(self.coboundaries),
            marked_loops = cp.deepcopy(self.marked_loops))
        if hasattr(self, "random"):
            Y.random = True
        return Y

    def orient(self):
        if self.oriented:
            return
        
        dim = self.dimension()
        bnd = self.boundaries      # a nested list three layers deep
        
        # in dimension 0 there are no boundaries
        # in dimension 1 orient them like [1, -1]
        for cell in bnd[1]:
            k = list(cell)
            assert(len(k) == 2)
            cell[k[0]] = 1
            cell[k[1]] = -1
        
        for d in range(2,dim+1):
            for cell in bnd[d]:
                if not cell:
                    continue
                for c in cell:
                    break

                cell[c] = 1                
                S = {c}
                T = set(cell)
                T.remove(c)

                while T:
                    cond = False;
                    add_to_S = []
                    for s in S:
                        for t in T:
                            L = set(bnd[d-1][s]) & set(bnd[d-1][t])
                            if not L:
                                continue
                            for l in L:
                                break
                            add_to_S.append(s)
                            T.remove(t)   
                            cond = True;
                            cell[t] = -cell[s] * bnd[d-1][s][l] * bnd[d-1][t][l]
                            break
                        if cond:
                            break
                    S.update(add_to_S)
        self.oriented = True

    def fix_bnd_loops(self):
        bnd = self.boundaries
        #rint([id(dc)-id(self._NULL_DICT) for dc in bnd[0]])
        perm = []  # perm[d][i] will be the new name of d-cell i

        for d in range(len(bnd)):
            perm.append(dict())
            skip = 0
            for n in range(len(bnd[d])):
                if bnd[d][n] is self._NULL_DICT:
                    skip += 1
                else:
                    perm[d][n] = n - skip
            #print(perm[d])
        
        for d in range(len(bnd)):
            bnd[d] = [x for x in bnd[d] if x is not self._NULL_DICT]
            if d > 0:
                for x in range(len(bnd[d])):
                    #try:
                    bnd[d][x] = {perm[d-1][c]: o for c,o in bnd[d][x].items()}
                    #except KeyError:
                    #    print(perm[0])
                        
        for loop in self.marked_loops:
            for i in range(len(loop)):
                e, s, t = loop[i]
                loop[i] = perm[1][e], perm[0][s], perm[0][t]
        
        return perm

    def cell_boundary(self, n, k):
        """
        Returns a set of pairs (d, i) indicating that d-cell i
        is in the boundary of n-cell k.
        """
        if n == 0:
            return set()

        N = n
        V = set(self.boundaries[N][k])   # V : set[int]
        cells = {(N-1, i) for i in V}
        N -= 1

        while N > 0:
            tmp = set()
            for v in V:
                tmp.update(self.boundaries[N][v])
                cells.update((N-1, i) for i in self.boundaries[N][v])
                #print(cells)
            V = tmp
            N -= 1
            
        return cells
    
    def join_cells(self, d, n):                
        #The n-th d-cell is removed assuming its coboundary has size 2
        bnd = self.boundaries     # list[list[dict[int, int]]]
        cobnd = self.coboundaries # list[list[set[int]]]
        # check the condition
        if not (len(cobnd[d][n]) == 2 and len(bnd[d][n]) > 0):
            return False
        c1, c2 = tuple(cobnd[d][n])
        if cobnd[d+1][c1] != cobnd[d+1][c2]:
            return False
        V1 = self.cell_boundary(d, n);
        V2 = self.cell_boundary(d+1, c1);
        V3 = self.cell_boundary(d+1, c2);
        if 1 + len(V1) != len(V2 & V3):  
            return False
        
        if d == 0:
            # we are merging two edges
            for loop in self.marked_loops:
                for i in range(len(loop)):
                    if loop[i] is None or loop[i-1] is None:
                        continue
                    if loop[i][0] == loop[i-1][0] == c1 or loop[i][0] == loop[i-1][0] == c2:
                        # =====> or <===== pattern
                        loop[i] = None
                        loop[i-1] = None
                    if loop[i][0] == c1 and loop[i-1][0] == c2:
                        _, s, t = loop[i-1]
                        loop[i] = c1, s, loop[i][2]
                        loop[i-1] = None
                        break
                    elif loop[i][0] == c2 and loop[i-1][0] == c1:
                        _, s, t = loop[i]
                        loop[i-1] = c1, loop[i-1][1], t
                        loop[i] = None
                        break
                loop2 = [x for x in loop if x is not None]
                del loop[:]
                loop.extend(loop2)
        elif d == 1 and self.marked_loops:
            # V2 and V3 contain only 1- and 0-cells so we can compare sizes
            V = V2 if len(V2) < len(V3) else V3
            V.remove((1, n)) # <---- the cell
            d01 = {c0: set() for dim, c0 in V if dim == 0}
            for dim, c1 in V:
                if dim != 1:
                    continue
                for c0 in bnd[1][c1]:
                    d01[c0].add(c1)
            k, h = tuple(bnd[1][n])
            kk = k
            # arbitrarily build path from s
            path = []
            while kk != h:
                for c in d01[kk]:
                    break
                for t in bnd[1][c]:
                    if t != kk:
                        break
                path.append((c, kk, t))
                d01[t].remove(c)
                d01[kk].remove(c)
                kk = t
            
            pathrev = [(e, t, s) for e, s, t in reversed(path)]
            
            for loop in self.marked_loops:
                loop2 = []
                for i in range(len(loop)):
                    if loop[i][0] == n:
                        if loop[i][1] == k:
                            loop2.extend(path)
                        else:
                            loop2.extend(pathrev)
                        break
                    else:
                        loop2.append(loop[i])
                del loop[:]
                loop.extend(loop2)
            
        # remove n from its boundaries
        for m in bnd[d][n]:
            cobnd[d-1][m].remove(n)

        # merge c2 into c1
        bnd[d][n] = self._NULL_DICT
        cobnd[d][n] = set()

        a = bnd[d+1][c1].pop(n)
        b = bnd[d+1][c2].pop(n)

        if self.oriented:
            for x in bnd[d+1][c2]:
                bnd[d+1][c2] *= -a*b
        
        bnd[d+1][c1].update(bnd[d+1][c2])

        # remove c2 from its coboundaries
        for m in cobnd[d+1][c2]:
            bnd[d+2][m].pop(c2)
        
        # remove c2 from its boundaries
        for m in bnd[d+1][c2]:
            if cobnd[d][m]:
                cobnd[d][m].remove(c2)
                cobnd[d][m].add(c1)
                
        bnd[d+1][c2] = self._NULL_DICT
        cobnd[d+1][c2] = set()
        return True    
    
    def simplify_once(self):
        bnd = self.boundaries     # bnd   : list[list[dict[int, int]]]
        cobnd = self.coboundaries # cobnd : list[list[set[int]]]

        for d in range(self.dimension()): # cant merge top dim cells
            for n in range(len(bnd[d])):
                if len(cobnd[d][n]) == 2 and len(bnd[d][n]) > 0:
                    self.join_cells(d,n)

        # cells joined, now change ALL names
        self.fix_bnd_loops()
        self.calculate_coboundaries()
    
    def simplify(self):
        a = self.size()
        b = -1
        while a > b:
            self.simplify_once()
            if b != -1:
                a = b
            b = self.size()
    
    def contract(self, n):
        #This function removes pairs of n- and (n+1)-cells if possible.

        if not hasattr(self, "vectorField") or self.vectorField == None:
            self.vectorField = [dict() for _ in range(self.dimension())]
            self.inverseVectorField = [dict() for _ in range(self.dimension())]
            self.bnd = cp.deepcopy(self.boundaries)
            self.cobnd = cp.deepcopy(self.coboundaries)

        if not hasattr(self, "FREE"):
            self.FREE = list(range(len(self.cobnd[n])))

        Free = [i for i in self.FREE if len(self.cobnd[n][i]) == 1]
        
        if hasattr(self, "random"):
            import random
            random.shuffle(Free)
        
        if not Free:
            del self.FREE
            return False

        idx = -1
        while True:
            idx += 1
            try:
                i = Free[idx]
            except IndexError:
                break
            if len(self.cobnd[n][i]) == 1:
                # remove cells (n, i) and (n+1, U)
                # where i is only in the boundary of U
                for U in self.cobnd[n][i]:
                    break
                self.vectorField[n][U] = i
                self.inverseVectorField[n][i] = U
                
                if n > 0:
                    b = self.bnd[n][i]
                    for j in b:
                        self.cobnd[n-1][j].remove(i)
                
                b = self.bnd[n+1][U]
                for j in b:
                    self.cobnd[n][j].remove(U)
                    if len(self.cobnd[n][j]) == 1:
                        Free.append(j)
                    
                  ###
                self.bnd[n][i] = self._NULL_DICT
                self.bnd[n+1][U] = self._NULL_DICT
                self.cobnd[n][i] = set()
                self.cobnd[n+1][U] = set()
                
                #print("Removing cells (%s, %s) and (%s, %s)"%(n,i,n+1,U))
                
                if n == 0:
                    # the only case here is if a loop retraces itself
                    for loop in self.marked_loops:
                        for j in range(len(loop)):
                            if loop[j] is not None and loop[j][1] == i:
                                loop[j] = None
                                loop[j-1] = None
                        loop2 = [ed for ed in loop if ed is not None]
                        del loop[:]
                        loop.extend(loop2)
                        
                elif n == 1:
                    # we need to replace every occurrence of this edge
                    # with the rest of the boundary of the boundary of U
                    V = self.cell_boundary(n+1, U);
                    V.remove((1, i)) # <---- the cell
                    
                    d01 = {c0: set() for dim, c0 in V if dim == 0}
                    #print(d01)
                    #print(V)
                    for dim, c1 in V:
                        if dim != 1:
                            continue
                        for c0 in self.boundaries[1][c1]:
                            #print("c1 = %s, d(c1) = %s, c0 = %s" % (c1, self.boundaries[1][c1], c0)) 
                            d01[c0].add(c1)
                    k, h = tuple(self.boundaries[1][i])
                    kk = k
                    # arbitrarily build path from s
                    path = []
                    while kk != h:
                        for c in d01[kk]:
                            break
                        for t in self.boundaries[1][c]:
                            if t != kk:
                                break
                        path.append((c, kk, t))
                        d01[t].remove(c)
                        d01[kk].remove(c)
                        kk = t
                    
                    #print(path)
                    pathrev = [(e, t, s) for e, s, t in reversed(path)]
                    
                    for loop in self.marked_loops:
                        loop2 = []
                        for i in range(len(loop)):
                            if loop[i][0] == n:
                                if loop[i][1] == k:
                                    loop2.extend(path)
                                else:
                                    loop2.extend(pathrev)
                                break
                            else:
                                loop2.append(loop[i])
                        del loop[:]
                        loop.extend(loop2)
                else:
                         # further deletions do not affect the 1-skeleton
                    pass # ^^^ did a geneticist write this?
        
        if Free:
            self.FREE = Free
            return True
        else:
            del self.FREE
            return False

    def contract_exhaust(self, d = 0):
        dim = self.dimension()
        contract_once = True
        was_contracted = True
        n = dim-1

        while was_contracted or n > d: 
            # this fails when n == d and it hasn't been contracted once
            was_contracted = False
            for n in reversed(range(d,dim)): # from dim-1 to d
                while contract_once:
                    contract_once = self.contract(n)
                    if contract_once:
                        was_contracted = True
                contract_once = True

    def contracted(self, dd = 0):
        Y = self.copy()
        Y.contract_exhaust(dd)

        Y.boundaries = Y.bnd
        del Y.bnd
        
        Y.perm = Y.fix_bnd_loops()
        Y.coboundaries = Y.calculate_coboundaries()
        return Y

    def delete_vertices(self, verts):
        rem = []
        rem.append(set(verts))
        for i in range(self.dimension()):
            cb = set()
            for j in rem[i]:
                cb.update(self.coboundaries[i][j])
            rem.append(cb)
        for i in range(self.dimension()+1):
            for j in rem[i]:
                self.boundaries[i][j] = self._NULL_DICT
        self.fix_bnd_loops()
        self.coboundaries = self.calculate_coboundaries()

    def cell_vertices(self, n, k):
        return {c for dim, c in self.cell_boundary(n,k) if dim == 0}

    def new_cell(self, n):
        self.boundaries[n].append(dict())
        return n, len(self.boundaries[n]) - 1

    def subdivide_edge(self, e):
        """
        Assumes that self is a triangulation.
        Destroys orientation data.
        """
        raise NotImplementedError()
        self.oriented = False
        # we need to add intermediate cells and double all affected cells.
        affected_cells = [None, {e}]
        for i in range(1, self.dimension()):
            affected_cells.append(set())
            for j in affected_cells[i]:
                affected_cells[i+1].update(self.coboundaries[i][j])
        
        a, b = tuple(self.boundaries[1][e])
        return
    
    def graded_coboundaries(self, n, k):
        gc = dict()
        gc[n] = {k}
        for d in range(n+1,self.dimension()+1):
            gc[d] = set()
            for j in gc[d-1]:
                gc[d].update(self.coboundaries[d-1][j])
                
        return gc
            
    def contract_edge(self, e):
        # condition: for every n-cell A in the coboundaries of a
        #      intersecting an n-cell B in the coboundaries of b
        # there exists an (n+1)-cell in the coboundaries of e
        # having A, B in its boundary
        raise NotImplementedError()
        a, b = tuple(self.boundaries[1][e])
        ca = self.graded_coboundaries(0,a)
        cb = self.graded_coboundaries(0,b)
        ce = self.graded_coboundaries(1,e)
        for d in range(1,self.dimension()+1):
            for A in ca[d]:
                for B in cb[d]:
                    if A == B:
                        continue
                    if self.boundaries[d][A] & self.boundaries[d][B]:
                        if self.coboundaries[d][A] & self.boundaries[d][B]:
                            pass
    
    def truncate(self, n):
        """Deletes all cells of dimension > n, effectively causing self.dimension() to return n."""
        del self.boundaries[n+1:]
        del self.coboundaries[n+1:]
        for s in self.coboundaries[n]:
            s.clear()
    
    def barycentric_subdivision(self):
        """Returns a list of tuples defining a simplicial complex homeomorphic to self."""
        _id_start = [0]
        for d in range(1,self.dimension()+1):
            _id_start.append(_id_start[-1] + self.ncells(d-1))
        def cell_id(n, k):
            """Returns an integer uniquely identifying a cell in self."""
            return _id_start[n] + k
        
        def enumerate_flag(n, k):
            """Returns an iterator of lists of cell_id values
            representing simplices in the barycentric subdivision of cell (n, k)."""
            if n == 0:
                yield [k]
                return
            else:
                idk = cell_id(n, k)
                for c in self.boundaries[n][k]:
                    for fl in enumerate_flag(n-1, c):
                        fl.append(idk)
                        yield fl
        
        return [tuple(sorted(fl)) for k in range(self.ncells(self.dimension())) for fl in enumerate_flag(self.dimension(), k)]
    
    def face_lattice_orientation(self):
        FL = [[] for d in range(self.dimension()+2)]
        ori = [[] for d in range(self.dimension()+2)]
        
        FL[0].extend([1,0] for _ in range(self.ncells(0)))
        ori[0].extend([1] for _ in range(self.ncells(0)))
        
        for d in range(1, self.dimension()+1):
            for k in range(len(self.boundaries[d])):
                ldk = list(self.boundaries[d][k].items())
                FL[d].append([len(ldk)] + [c+1 for c,o in ldk])
                ori[d].append([o for c,o in ldk])
        
        return FL, ori
