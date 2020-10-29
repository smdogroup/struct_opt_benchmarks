import numpy as np
import openmdao.api as om
import pickle
import argparse
from paropt.paropt_driver import ParOptDriver
from plane_stress_analysis import PlaneStressAnalysis
import os
import timeit

# Set up parameters that we want to vary
p = argparse.ArgumentParser()
p.add_argument('picklename', metavar='[picklename]', type=str)
p.add_argument('--optimizer', default='ParOpt',
                    choices=['ParOpt', 'SNOPT', 'IPOPT'])
p.add_argument('--freq', type=float, default=0.5,
                    help='lower bound of natural frequency')
p.add_argument('--mass', type=float, default=0.4,
                    help='upper bound of mass, in fraction of full mass')
p.add_argument('--max_iter', type=int, default=200,
                    help='maximum number of major iteration')
p.add_argument('--ParOpt_use_adaptive_gamma_update', action='store_true')
p.add_argument('--ParOpt_use_filter', action='store_true')
p.add_argument('--ParOpt_use_soc', action='store_true')
p.add_argument('--info', type=str, default='')
p.add_argument('--outdir', type=str, default='')
args = p.parse_args()

# create directory if outdir doesn't exist
if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)

# Set up constants that we want to fix
eigshsigma = -100.0  # sigma for eigsh's shift-invert mode
eig_ks = 100.0  # ks parameter for minimum eigenvalue calculation
num_eigs = 8  # number of eigenvalues computed
qn_diag_type = 'yts_over_sts'
x_init = 0.95  # initial guess for design variable

# Save optimization-specific settings for postprocess' purpose
opt_settings = dict()
opt_settings['eigshsigma'] = eigshsigma
opt_settings['eig_ks'] = eig_ks
opt_settings['num_eigs'] = num_eigs
opt_settings['freq'] = args.freq

# Load in pickle file
with open(args.picklename, 'rb') as pklfile:
    prob_pkl = pickle.load(pklfile)

# Get preprocessed problem
prob_name = prob_pkl['prob_name']
nelems = prob_pkl['nelems']
nnodes = prob_pkl['nnodes']
ndof = prob_pkl['ndof']
C = prob_pkl['C']
conn = prob_pkl['conn']
X = prob_pkl['X']
dof = prob_pkl['dof']
force = prob_pkl['force']
r0 = prob_pkl['r0']
density = prob_pkl['density']
qval = prob_pkl['qval']
x = prob_pkl['x']

# Initialize analysis object
lambda0 = (2.0*np.pi*float(args.freq))**2
analysis = PlaneStressAnalysis(conn, dof, X, force,
    r0, qval, C, density=density, freqconstr=True,
    lambda0=lambda0, ks_parameter=eig_ks, num_eigs=num_eigs,
    eigshsigma=eigshsigma)

# Compute mass for a fully-filled structure
xfull = np.ones(nnodes)
full_mass = analysis.mass(xfull)

# Form openMDAO problem
prob = om.Problem()
x = x_init*np.ones(nnodes)
# np.random.seed(0)
# x = np.random.rand(nnodes)
indeps = prob.model.add_subsystem('indeps', om.IndepVarComp())
indeps.add_output('x', x)
prob.model.add_subsystem('topo', analysis)
prob.model.connect('indeps.x', 'topo.x')
prob.model.add_design_var('indeps.x', lower=1e-3, upper=1.0)
prob.model.add_objective('topo.c', scaler=1.0)
prob.model.add_constraint('topo.m', upper=full_mass*args.mass)
prob.model.add_constraint('topo.freq', lower=0.0)

# Setup output format
extra = ''
if args.optimizer == 'ParOpt':
    if args.ParOpt_use_adaptive_gamma_update:
        extra += 'adapt'
    if args.ParOpt_use_filter:
        extra += 'filter'
    if args.ParOpt_use_soc:
        extra += 'soc'
info = ''
info += args.info
if info != '':
    info = '-' + info
outputname = 'comp_min_massfreq_constr-{:s}-mass-{:.3f}-freq-{:.3f}-{:s}{:s}'.format(
    args.optimizer+extra, args.mass, args.freq, prob_name, info)
if args.outdir != '':
    outputname = args.outdir + '/' + outputname

# Setup optimizer
if args.optimizer == 'ParOpt':
    if args.ParOpt_use_filter:
        tr_accept_step_strategy = 'filter_method'
    else:
        tr_accept_step_strategy = 'penalty_method'
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
        'tr_l1_tol': 0.0,
        'tr_linfty_tol': 0.0,
        'tr_adaptive_gamma_update':
            args.ParOpt_use_adaptive_gamma_update,
        'tr_accept_step_strategy':
            tr_accept_step_strategy,
        'filter_sufficient_reduction': True,
        'tr_use_soc': args.ParOpt_use_soc,
        'tr_max_iterations': args.max_iter,
        'output_file': outputname+'.out',
        'tr_output_file': outputname+'.tr',
        'penalty_gamma': 50.0,
        'qn_subspace_size': 2,
        'qn_type': 'bfgs',
        'qn_diag_type': qn_diag_type,
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
    # prob.driver.opt_settings['file_print_level'] = 6

elif args.optimizer == 'SNOPT':
    prob.driver = om.pyOptSparseDriver()
    prob.driver.options['optimizer'] = 'SNOPT'
    prob.driver.opt_settings['Iterations limit'] = 9999999999999
    prob.driver.opt_settings['Major iterations limit'] = args.max_iter
    prob.driver.opt_settings['Summary file'] = outputname+'.out'
    prob.driver.opt_settings['Print file'] = outputname+'_print.out'
    prob.driver.opt_settings['Major print level'] = 1
    prob.driver.opt_settings['Minor print level'] = 0

#  Run optimization
t_start = timeit.default_timer()
prob.setup()
prob.run_driver()
t_end = timeit.default_timer()

# Get design
x_opt = prob.get_val('indeps.x')
analysis.plot_solution(x_opt, savefig=True, name=outputname)

# save result to solution pickle file
prob_pkl['x'] = x_opt
prob_pkl['opt_settings'] = opt_settings
prob_pkl['opt_time'] = '{:.2e} s'.format(t_end - t_start)
picklename = outputname+'.pkl'
with open(picklename, 'wb') as pklfile:
    pickle.dump(prob_pkl, pklfile)
