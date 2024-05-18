"""
Contains a Python adaptation and extension of Regina source code.

For performance reasons, barycentric subdivisions of 5-dimensional (and higher)
triangulations are disabled by default in Regina.
This is a Python rewriting of the same code without the restriction.

Regina is released under GPL 2.0.
Copyright (c) 1999–2023, The Regina development team.
"""
from regina import *

__all__ = ["bary_sub5"]

def bary_sub5(tri):
    """
    Computes the barycentric subdivision of a Triangulation5 object.
    
    Exceptions:
    
        - ValueError: if "tri" is not a Triangulation5.
    """
    
    if not isinstance(tri, Triangulation5):
        raise ValueError("tri must be a Triangulation5 instance")
        
    nOld = tri.fVector()[-1]
    if nOld == 0:
        return tri

    staging = Triangulation5()
    PERMS_6 = 720
    
    # A top-dimensional simplex in the subdivision is uniquely defined
    # by a permutation p on (5+1) elements.
    #
    # As described in the documentation for barycentricSubdivision(),
    # this is the simplex that:
    # - meets the boundary in the facet opposite vertex p[5];
    # - meets that facet in the (5-2)-face opposite vertex p[5-1];
    # - meets that (5-2)-face in the (5-3)-face opposite vertex p[5-2];
    # - ...
    # - meets that edge in the vertex opposite vertex p[1];
    # - directly touches vertex p[0].

    newSimp = [staging.newSimplex() for _ in range(PERMS_6*nOld)]

    # Do all of the internal gluings
    permIdx = adjIdx = 0
    perm = Perm6()
    glue = Perm6()
    i = 0
    
    vdict = {} # old vertex -> new simplex
    
    for simp in range(nOld):
        for permIdx in range(PERMS_6):
            perm = Perm6.orderedSn[permIdx]
            
            # Internal gluings within the old simplex:
            for i in range(5):
                adjIdx = (perm * Perm6(i, i+1)).orderedSnIndex()
                if (permIdx < adjIdx):
                    newSimp[PERMS_6 * simp + permIdx].join(perm[i],
                        newSimp[PERMS_6 * simp + adjIdx],
                        Perm6(perm[i], perm[i+1]))
            # the zeroth vertices of the new simplices
            # are the original vertices

            # Adjacent gluings to the adjacent simplex:
            oldSimp = tri.simplex(simp)
            if not oldSimp.adjacentSimplex(perm[5]):
                continue # This hits a boundary facet.
            if newSimp[PERMS_6 * simp + permIdx].adjacentSimplex(perm[5]):
                continue # We've already done this gluing from the other side.
            glue = oldSimp.adjacentGluing(perm[5])
            newSimp[PERMS_6 * simp + permIdx].join(perm[5],
                newSimp[PERMS_6 * oldSimp.adjacentSimplex(perm[5]).index() + (glue * perm).orderedSnIndex()],
                glue)
    
    # vertex I 
    
    return staging
