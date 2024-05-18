"""
Contains a Python adaptation and extension of Regina source code.

For performance reasons, barycentric subdivisions of 5-dimensional (and higher)
triangulations are disabled by default in Regina.
This is a Python rewriting of the same code without the restriction.

Regina is released under GPL 2.0.
Copyright (c) 1999–2023, The Regina development team.
"""
import regina
from sage.all import factorial

__all__ = ["bary_sub"]

def bary_sub(tri):
    """
    Computes the barycentric subdivision of a Triangulation object.
    
    Exceptions:
    
        - ValueError: if "tri" is not a Triangulation.
    """
    DIM = tri.dimension
    TriClass = getattr(regina, f"Triangulation{DIM}")
    PermClass = getattr(regina, f"Perm{DIM+1}")
    if not isinstance(tri, TriClass):
        raise ValueError(f"tri must be a Triangulation{DIM} instance")
        
    nOld = tri.fVector()[-1]
    if nOld == 0:
        return tri

    staging = TriClass()
    PERMS = int(factorial(DIM+1))
    
    # A top-dimensional simplex in the subdivision is uniquely defined
    # by a permutation p on (dim+1) elements.
    #
    # As described in the documentation for barycentricSubdivision(),
    # this is the simplex that:
    # - meets the boundary in the facet opposite vertex p[dim];
    # - meets that facet in the (dim-2)-face opposite vertex p[dim-1];
    # - meets that (dim-2)-face in the (dim-3)-face opposite vertex p[dim-2];
    # - ...
    # - meets that edge in the vertex opposite vertex p[1];
    # - directly touches vertex p[0].

    newSimp = [staging.newSimplex() for _ in range(PERMS*nOld)]

    # Do all of the internal gluings
    permIdx = adjIdx = 0
    perm = PermClass()
    glue = PermClass()
    i = 0
    
    vdict = {} # old vertex -> new simplex
    
    for simp in range(nOld):
        for permIdx in range(PERMS):
            perm = PermClass.orderedSn[permIdx]
            
            # Internal gluings within the old simplex:
            for i in range(DIM):
                adjIdx = (perm * PermClass(i, i+1)).orderedSnIndex()
                if (permIdx < adjIdx):
                    newSimp[PERMS * simp + permIdx].join(perm[i],
                        newSimp[PERMS * simp + adjIdx],
                        PermClass(perm[i], perm[i+1]))
            # the zeroth vertices of the new simplices
            # are the original vertices

            # Adjacent gluings to the adjacent simplex:
            oldSimp = tri.simplex(simp)
            if not oldSimp.adjacentSimplex(perm[DIM]):
                continue # This hits a boundary facet.
            if newSimp[PERMS * simp + permIdx].adjacentSimplex(perm[DIM]):
                continue # We've already done this gluing from the other side.
            glue = oldSimp.adjacentGluing(perm[DIM])
            newSimp[PERMS * simp + permIdx].join(perm[DIM],
                newSimp[PERMS * oldSimp.adjacentSimplex(perm[DIM]).index() + (glue * perm).orderedSnIndex()],
                glue)
    
    # vertex I 
    
    return staging
