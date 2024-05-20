class FGARankOptions:
    """
    Configuration for the computations of ranks over the algebra of a finite group.
    
        - BICYCLIC:
            When reducing a finite group algebra matrix, invert bicyclic units.
            Default: True.
            
        - TRANSPOSE:
            Transpose before calculating rank.
            Default: False.
            
        - CEILING:
            Return the ceiling of the rank.
            Default: False.
            
        - SPARSE:
            Use sparse matrices.
            Default: True.
            
        - RING:
            Ring used for rank computations. Allowed values:

            - "rational" : rationals (fastest, exact).
            - "real"     : reals with 53 bits of precision, inaccurate due to rounding errors.
            - "double"   : reals with NumPy, probably accurate but slow, at least for small matrices.
            
            Default: "rational".
            
        - TRY_REDUCE:
            Use Gauss moves before restriction of scalars.
            Default: False.
    """
    
    BICYCLIC = True
    TRANSPOSE = False
    CEILING = False
    SPARSE = True
    RING = "rational"
    TRY_REDUCE = False
    
class TwistedRankOptions:
    """
    Configuration for the von_neumann_rank and determinant_degree routines.
    
        - PRINT_FIN:
            Print the matrix over Q[L].
            Default: False.
            
        - ASK_INPUT:
            Ask whether to save a copy of the expanded matrix and the matrix over Q[L].
            Default: False.
    """
    
    PRINT_FIN = False
    ASK_INPUT = False

class LogOptions:
    """
    Configuration for logging.
    
        - LEVEL:
            Level of logging verbosity.
            Default: SILENT.
        
        - PRECISION:
            Precision of values displayed during logging (number of digits, or PREC_EXACT).
            Default: 15.
            
    Constants:
        
        - PREC_EXACT:
            Constant value for `PRECISION`. Print rational values as fractions.
    """
    
    LEVEL = 0
    PREC_EXACT = -1
    PRECISION = 15
