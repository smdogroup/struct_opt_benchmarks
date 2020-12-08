# First, import fortran analysis libraries
try:
    import plane_lib
except ImportError:
    raise ImportError("\nPlane stress analysis library is missing, please compile fortran code first!\n")

try:
    import solid_lib
except ImportError:
    raise ImportError("\nSolid analysis library is missing, please compile fortran code first!\n")

# Next, import other modules
from py_struct_opt_benchmarks import wrapper
from py_struct_opt_benchmarks import utils
from py_struct_opt_benchmarks import preprocessors
from py_struct_opt_benchmarks import optproblems

