# First, import fortran analysis library
try:
    from plane_stress import plane_lib
except ImportError:
    raise ImportError("\nAnalysis library is missing, please compile fortran code first!\n")

# Next, import other modules
from plane_stress import wrapper
from plane_stress import utils
from plane_stress import preprocessors
from plane_stress import optproblems

