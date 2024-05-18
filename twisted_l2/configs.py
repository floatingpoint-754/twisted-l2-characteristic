class FGARankOptions:
    """When reducing a finite group algebra matrix, invert bicyclic units."""
    BICYCLIC = True
    
    """Transpose before calculating rank."""
    TRANSPOSE = False
    
    """Return the ceiling of the rank."""
    CEILING = False
    
    """Use sparse matrices."""
    SPARSE = True
    
    """
    Ring used for rank computations. Allowed values:
    
        - "rational" : rationals (fastest, exact).
        - "real"     : reals with 53 bits of precision, inaccurate due to rounding errors.
        - "double"   : reals with NumPy, probably accurate but slow, at least for small matrices.
    """
    RING = "rational"
    
    """Use Gauss moves before restriction of scalars."""
    TRY_REDUCE = False
    
class TwistedRankOptions:
    """Print the matrix over Q[L]."""
    PRINT_FIN = False
    
    """Ask whether to save a copy of the expanded matrix and the matrix over Q[L]."""
    ASK_INPUT = False

class LogOptions:
    """Level of logging verbosity"""
    LEVEL = 0
    
    """Constant value for `PRECISION`. Print rational values as fractions."""
    PREC_EXACT = -1
    
    """Precision of values displayed during logging (number of digits, or `PREC_EXACT`)"""
    PRECISION = 15
