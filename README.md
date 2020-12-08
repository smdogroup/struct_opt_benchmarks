In order to use this repo as a python package, you need to add /struct_opt_benchmarks in $PYTHONPATH

The work flow using this repo to solve topology optimization problems will be:

1. Use `preprocessors/preproc.py` to generate a .pkl file that contains the description of entire problem,
including geometry, mesh, load and boundary condition. Run `preproc.py --help` to see all options.
Note that tmr needs to be installed in order to generate unstructured mesh. Please refer to tmr document for
more information: [https://smdogroup.github.io/tmr/](https://smdogroup.github.io/tmr/).

2. Execute `optproblems/optimize.py` with the problem pkl file to run optimizations. Run `optimize.py --help`
to see all options, including all supported optimization problems, optimizer settings, etc.

3. Use utility functions in `utils/` to postprocess results.
