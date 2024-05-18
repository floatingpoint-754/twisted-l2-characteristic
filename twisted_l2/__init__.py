"""
This package provides methods for the computation of the twisted L^2-Euler characteristic.

The algorithm roughly works as follows:

1. start with a finite connected L^2-acyclic CW complex Y,
   with residually finite fundamental group G satisfying the Atiyah conjecture;
  
2. construct the cellular chain complex of the universal cover of Y
   (which can be seen as a chain complex over Z[G]);
  
3. given a homomorphism phi: G --> Z, the twisted L^2-Euler characteristic is
   a certain sum involving the combinatorial Laplacians of the chain complex;
  
4. more precisely, we take the "phi-degree" of the Dieudonné determinant
   of each Laplacian;
  
5. this "phi-degree" is related to the order valuation, and can be calculated
   using Oki's matrix expansion algorithm;
  
6. this reduces the problem to computing von Neumann ranks of matrices,
   which can be approximated (by Lück's theorem) by ranks of certain rational matrices.
  
Step 6 is usually the most computationally expensive, due to the sheer size
of the matrices involved. However, for some very complicated spaces Y
(e.g. the Ratcliffe-Tschantz manifold), the real bottleneck is step 2,
which does not even run to completion in reasonable time.

This package depends on the packages "sage" and "regina".
"""

from .twisted_l2_char import *
from .twisted_rank import *
from .utils import *
from .serializer import *
from .free_by_cyclic import *
from .fga_rank import *
from .equivariant_boundary import *
from .logger import *

