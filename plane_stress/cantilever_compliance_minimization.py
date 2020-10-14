from mpi4py import MPI
import numpy as np
import argparse
import openmdao.api as om
from paropt.paropt_driver import ParOptDriver
import re

from plane_stress_analysis import PlaneStressAnalysis

# Define command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--optimizer', default='ParOpt',
                    choices=['ParOpt', 'SNOPT', 'IPOPT'])
parser.add_argument('--freq_constr', action='store_true',
                    help='switch on frequency constraint')
parser.add_argument('--freq', type=float, default=0.5,
                    help='minimum frequency required')
parser.add_argument('--max_iter', type=int, default=200,
                    help='maximum major iteration')
parser.add_argument('--x_init', type=str, default=None,
                    help='initial design, should be a .npy file')
args = parser.parse_args()

# Define design space geometry, meshing, material property
# and boundary conditions
nx = 64
ny = 64
r0 = 1.0/32.0
lx = 1.0
ly = 1.0
forceval = 25.0

nelems = nx*ny
nnodes = (nx+1)*(ny+1)

conn = np.zeros((nelems, 4), dtype=np.intc)
vars = -np.ones((nnodes, 2), dtype=np.intc)
X = np.zeros((nnodes, 2))
rho = np.ones(nnodes)

C = np.zeros((3, 3))
qval = 5.0

density = 2700.0 # kg/m^3

frequency = args.freq
lambda0 = (2.0*np.pi*frequency)**2
ks_parameter = 100.0

E = 70e3
nu = 0.3
C[0, 0] = E/(1.0 - nu**2)
C[0, 1] = nu*E/(1.0 - nu**2)
C[1, 0] = C[0, 1]
C[1, 1] = C[0, 0]
C[2, 2] = 0.5*E/(1.0 + nu)

for j in range(ny):
    for i in range(nx):
        conn[i + j*nx, 0] = i + (nx+1)*j
        conn[i + j*nx, 1] = i+1 + (nx+1)*j
        conn[i + j*nx, 2] = i + (nx+1)*(j+1)
        conn[i + j*nx, 3] = i+1 + (nx+1)*(j+1)

nvars = 0
for j in range(ny+1):
    for i in range(nx+1):
        X[i + j*(nx+1), 0] = lx*i/nx
        X[i + j*(nx+1), 1] = ly*j/ny
        if i > 0:
            vars[i + j*(nx+1), 0] = nvars
            nvars += 1
            vars[i + j*(nx+1), 1] = nvars
            nvars += 1

force = np.zeros(nvars)
i = nx
j = 0
force[vars[i + j*(nx+1), 1]] = -forceval

# Setup openmdao problem object
prob = om.Problem()

# Create independent variable component
x_init = 0.95*np.ones(nnodes)
indeps = prob.model.add_subsystem('indeps', om.IndepVarComp())
indeps.add_output('x', x_init)

# Add analysis
analysis = PlaneStressAnalysis(
    conn, vars, X, force, r0, qval, C, density,
    freqconstr=args.freq_constr, lambda0=lambda0)
prob.model.add_subsystem('topo', analysis)
prob.model.connect('indeps.x', 'topo.x')

# Form optimization problem
prob.model.add_design_var('indeps.x', lower=1e-3, upper=1.0)
prob.model.add_objective('topo.c', scaler=1.0)
prob.model.add_constraint('topo.m', lower=0.0)
if args.freq_constr is True:
    prob.model.add_constraint('topo.freq', lower=0.0)

# Set up output format
if args.freq_constr is True:
    outputname = '{:s}-freq-{:.1f}-iter-{:d}'.format(
        args.optimizer, args.freq, args.max_iter)
else:
    outputname = '{:s}-iter-{:d}'.format(
        args.optimizer, args.max_iter)

# Set up optimizer
if args.optimizer == 'ParOpt':
    prob.driver = ParOptDriver()
    options = {
        'algorithm': 'tr',
        'output_level':0,
        'norm_type': 'l1',
        'tr_init_size': 0.05,
        'tr_min_size': 0.001,
        'tr_max_size': 10.0,
        'tr_eta': 0.25,
        'tr_infeas_tol': 1e-6,
        'tr_l1_tol': 1e-3,
        'tr_linfty_tol': 0.0,
        'tr_adaptive_gamma_update': False,
        'tr_max_iterations': args.max_iter,
        'output_file': outputname+'.out',
        'tr_output_file': outputname+'.tr',
        'penalty_gamma': 50.0,
        'qn_subspace_size': 2,
        'qn_type': 'bfgs',
        'qn_diag_type': 'yts_over_sts',
        'abs_res_tol': 1e-8,
        'starting_point_strategy': 'affine_step',
        'barrier_strategy': 'mehrotra_predictor_corrector',
        'tr_steering_barrier_strategy':
            'mehrotra_predictor_corrector',
        'tr_steering_starting_point_strategy': 'affine_step',
        'use_line_search': False}
    for key in options:
        prob.driver.options[key] = options[key]

elif args.optimizer == 'IPOPT':
    prob.driver = om.pyOptSparseDriver()
    prob.driver.options['optimizer'] = 'IPOPT'
    prob.driver.opt_settings['max_iter'] = args.max_iter
    prob.driver.opt_settings['output_file'] = outputname+'.out'


elif args.optimizer == 'SNOPT':
    prob.driver = om.pyOptSparseDriver()
    prob.driver.options['optimizer'] = 'SNOPT'
    prob.driver.opt_settings['Iterations limit'] = 9999999999999
    prob.driver.opt_settings['Major iterations limit'] = args.max_iter
    prob.driver.opt_settings['Print file'] = outputname+'_print.out'
    prob.driver.opt_settings['Summary file'] = outputname+'_summary.out'
    prob.driver.opt_settings['Major print level'] = 1
    prob.driver.opt_settings['Minor print level'] = 0


# Run problem
prob.setup()
prob.run_driver()

# Output result
xopt = prob.get_val('indeps.x')
if args.freq_constr is True:
    np.save('{:s}-xopt-freq-{:.1f}-iter-{:d}'.format(
        args.optimizer, args.freq, args.max_iter), xopt)
    analysis.plot_solution(xopt, savefig=True,
        name='{:s}-freq-{:.1f}-iter-{:d}'.format(
        args.optimizer, args.freq, args.max_iter))
else:
    np.save('{:s}-xopt-iter-{:d}'.format(
        args.optimizer, args.max_iter), xopt)
    analysis.plot_solution(xopt, savefig=True,
        name='{:s}-iter-{:d}'.format(
        args.optimizer, args.max_iter))

