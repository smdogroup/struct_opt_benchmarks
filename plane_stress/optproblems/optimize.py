#!/usr/bin/env python3

"""
This script executes one of the following the topology optimization problems:
    - [comp_min_mass_constr]: compliance minimization, mass constrained
    - [comp_min_massfreq_constr]: compliance minimization, mass and frequency constrained
    - [comp_min_massstress_constr]: compliance minimization, mass and stress constrained
    - [comp_min_massfreqstress_constr]: compliance minimization, mass, stress and freq constrained
    - [stress_min_mass_constr] : stress minimization, mass constrained
    - [mass_min_stress_constr]: mass minimization, stress constrained
"""

import numpy as np
import openmdao.api as om
import pickle
import argparse
from paropt.paropt_driver import ParOptDriver
from plane_stress.wrapper import PlaneStressAnalysis
import os
import timeit

def optimize(in_picklename, in_opt_problem, in_optimizer, in_qvals, in_epsilon,
             in_ks_parameter, in_design_mass, in_design_freq, in_design_stress,
             in_stress_as_fraction, in_max_iters, in_ParOpt_use_adaptive_gamma_update,
             in_ParOpt_use_filter, in_ParOpt_use_soc, in_info, in_outdir):

    # Check inputs
    opt_choices = [
        'comp_min_mass_constr', 'comp_min_massfreq_constr',
        'comp_min_massstress_constr', 'comp_min_massfreqstress_constr',
        'stress_min_mass_constr', 'mass_min_stress_constr']

    if in_opt_problem not in opt_choices:
        raise ValueError("\nopt_problem {:s} is not supported!".format(in_opt_problem))

    if in_optimizer not in ['ParOpt', 'SNOPT', 'IPOPT']:
        raise ValueError("\noptimizer {:s} not supported!".format(in_optimizer))

    if len(in_qvals) != len(in_max_iters):
        raise ValueError("\nNumber of qvals does not match number of max_iter!")

    # Create directory if outdir doesn't exist
    if in_outdir is not None:
        # Check if specified output directory exists
        if os.path.isdir(in_outdir):
            # in_outdir exists
            pass
        # Otherwise, try to create the folder
        else:
            try:
                os.mkdir(in_outdir)
            except:
                print("\n[Warning] Cannot create directory {:s}!".format(in_outdir))
                in_outdir = None

    # Switch on continuation if more than one qval/max_iter specified
    if len(in_qvals) > 1 or len(in_max_iters) > 1:
        print("[Notice] more than one qval/max_iter is specified, using continuation...")

    # Set up constants
    x_init = 0.95
    ParOpt_qn_diag_type = 'yty_over_yts'
    ParOpt_penalty_gamma = 50.0
    qn_subspace_size = 5

    # Load pickle file
    with open(in_picklename, 'rb') as pklfile:
        prob_pkl = pickle.load(pklfile)

    # Get data from pickle
    prob_name = prob_pkl['prob_name']
    nnodes = prob_pkl['nnodes']
    C = prob_pkl['C']
    conn = prob_pkl['conn']
    X = prob_pkl['X']
    dof = prob_pkl['dof']
    force = prob_pkl['force']
    r0 = prob_pkl['r0']
    density = prob_pkl['density']

    # Set up info string for output
    optimizer_info = in_optimizer

    if in_optimizer == 'ParOpt':
        if in_ParOpt_use_adaptive_gamma_update:
            optimizer_info += 'Adapt'
        if in_ParOpt_use_filter:
            optimizer_info += 'Filter'
        if in_ParOpt_use_soc:
            optimizer_info += 'Soc'
    prob_info = ''

    # we might use continuation strategy to update qval and restart
    # optimization if more than one qval is specified, in which case
    # we use loop over all qvals
    index = 0
    x_temp = None
    for in_qval in in_qvals:

        print("[Notice] continuation step {:d}/{:d}".format(index+1, len(in_qvals)))

        # Get the number of iterations for current continuation step
        in_max_iter = in_max_iters[index]

        # Set up optimization problem
        prob = om.Problem()

        # If this is the first continuation step, use default starting point
        if index == 0:
            x = x_init*np.ones(nnodes)
        else:
            x = x_temp

        indeps = prob.model.add_subsystem('indeps', om.IndepVarComp())
        indeps.add_output('x', x)

        if in_opt_problem == 'comp_min_mass_constr':

            # Check inputs
            if in_design_mass is None:
                raise ValueError("\ndesign_mass must be specified!\n")

            # Create analysis object
            analysis = PlaneStressAnalysis(conn, dof, X, force, r0,
                C, density, in_qval, in_epsilon, in_ks_parameter,
                compute_comp=True, compute_mass=True)

            # Define openmdao problem
            prob.model.add_subsystem('topo', analysis)
            prob.model.connect('indeps.x', 'topo.x')
            x_full = np.ones(nnodes)
            full_mass = analysis.mass(x_full)
            prob.model.add_design_var('indeps.x', lower=1e-3, upper=1.0)
            prob.model.add_objective('topo.c', scaler=1.0)
            prob.model.add_constraint('topo.m', upper=full_mass*in_design_mass)

            # Set up problem info for output
            prob_info = 'mass-{:.2f}'.format(in_design_mass)

        elif in_opt_problem == 'comp_min_massfreq_constr':

            # Check inputs
            if in_design_mass is None:
                raise ValueError("\ndesign_mass must be specified!\n")

            if in_design_freq is None:
                raise ValueError("\ndesign_freq must be specified!\n")

            # Save to pkl file for postprocess' purpose
            prob_pkl['design_freq'] = in_design_freq

            # Create analysis object
            analysis = PlaneStressAnalysis(conn, dof, X, force, r0,
                C, density, in_qval, in_epsilon, in_ks_parameter,
                compute_comp=True, compute_mass=True, compute_freq=True,
                design_freq=in_design_freq)

            # Define openmdao problem
            prob.model.add_subsystem('topo', analysis)
            prob.model.connect('indeps.x', 'topo.x')
            x_full = np.ones(nnodes)
            full_mass = analysis.mass(x_full)
            prob.model.add_design_var('indeps.x', lower=1e-3, upper=1.0)
            prob.model.add_objective('topo.c', scaler=1.0)
            prob.model.add_constraint('topo.m', upper=full_mass*in_design_mass)
            prob.model.add_constraint('topo.freqc', lower=0.0)

            # Set up problem info for output
            prob_info = 'mass-{:.2f}-freq-{:.3f}'.format(in_design_mass, in_design_freq)

        elif in_opt_problem == 'comp_min_massstress_constr':

            # Check inputs
            if in_design_stress is None:
                raise ValueError("\ndesign_stress must be specified!\n")

            if in_design_mass is None:
                raise ValueError("\ndesign_mass must be specified!\n")

            if in_stress_as_fraction:

                analysis = PlaneStressAnalysis(conn, dof, X, force, r0,
                    C, density, in_qval, in_epsilon, in_ks_parameter)
                x_full = np.ones(nnodes)
                stress_full = np.max(analysis.nodal_stress(x_full))
                design_stress = in_design_stress*stress_full

            else:
                design_stress = in_design_stress

            # Save to pkl file for postprocess' purpose
            prob_pkl['design_stress'] = design_stress

            # Create analysis object
            analysis = PlaneStressAnalysis(conn, dof, X, force, r0,
                    C, density, in_qval, in_epsilon, in_ks_parameter,
                    compute_comp=True, compute_stress=True,
                    compute_mass=True, design_stress=design_stress)

            # Define openmdao problem
            prob.model.add_subsystem('topo', analysis)
            prob.model.connect('indeps.x', 'topo.x')
            x_full = np.ones(nnodes)
            full_mass = analysis.mass(x_full)
            prob.model.add_design_var('indeps.x', lower=1e-3, upper=1.0)
            prob.model.add_objective('topo.c', scaler=1.0)
            prob.model.add_constraint('topo.ks_nodal_stress', upper=1.0)
            prob.model.add_constraint('topo.m', upper=full_mass*in_design_mass)

            # Set up problem info for output
            prob_info = 'stress-{:.2e}-mass-{:.2f}'.format(design_stress, in_design_mass)

        elif in_opt_problem == 'comp_min_massfreqstress_constr':

            # Check inputs
            if in_design_stress is None:
                raise ValueError("\ndesign_stress must be specified!\n")

            if in_design_mass is None:
                raise ValueError("\ndesign_mass must be specified!\n")

            if in_design_freq is None:
                raise ValueError("\ndesign_freq must be specified!\n")

            if in_stress_as_fraction:

                analysis = PlaneStressAnalysis(conn, dof, X, force, r0,
                    C, density, in_qval, in_epsilon, in_ks_parameter)
                x_full = np.ones(nnodes)
                stress_full = np.max(analysis.nodal_stress(x_full))
                design_stress = in_design_stress*stress_full

            else:
                design_stress = in_design_stress

            # Save to pkl file for postprocess' purpose
            prob_pkl['design_stress'] = design_stress
            prob_pkl['design_freq'] = in_design_freq

            # Create analysis object
            analysis = PlaneStressAnalysis(conn, dof, X, force, r0,
                C, density, in_qval, in_epsilon, in_ks_parameter, compute_comp=True,
                compute_stress=True, compute_freq=True, compute_mass=True,
                design_stress=design_stress, design_freq=in_design_freq)

            # Define openmdao problem
            x_full = np.ones(nnodes)
            full_mass = analysis.mass(x_full)
            prob.model.add_subsystem('topo', analysis)
            prob.model.connect('indeps.x', 'topo.x')
            prob.model.add_design_var('indeps.x', lower=1e-3, upper=1.0)
            prob.model.add_objective('topo.c', scaler=1.0)
            prob.model.add_constraint('topo.ks_nodal_stress', upper=1.0)
            prob.model.add_constraint('topo.m', upper=full_mass*in_design_mass)
            prob.model.add_constraint('topo.freqc', lower=0.0)

            # Set up problem info for output
            prob_info = 'stress-{:.2e}'.format(design_stress)

        elif in_opt_problem == 'stress_min_mass_constr':

            # Check inputs
            if in_design_mass is None:
                raise ValueError("\ndesign_mass must be specified!\n")

            if in_design_stress is None:
                raise ValueError("\ndesign_stress must be specified!\n")

            if in_stress_as_fraction:

                analysis = PlaneStressAnalysis(conn, dof, X, force, r0,
                    C, density, in_qval, in_epsilon, in_ks_parameter)
                x_full = np.ones(nnodes)
                stress_full = np.max(analysis.nodal_stress(x_full))
                design_stress = in_design_stress*stress_full

            else:
                design_stress = in_design_stress

            # Save to pkl file for postprocess' purpose
            prob_pkl['design_stress'] = design_stress

            # Create analysis object
            analysis = PlaneStressAnalysis(conn, dof, X, force, r0,
                    C, density, in_qval, in_epsilon, in_ks_parameter,
                    compute_stress=True, compute_mass=True, design_stress=design_stress)

            # Define openmdao problem
            prob.model.add_subsystem('topo', analysis)
            prob.model.connect('indeps.x', 'topo.x')
            x_full = np.ones(nnodes)
            full_mass = analysis.mass(x_full)
            prob.model.add_design_var('indeps.x', lower=1e-3, upper=1.0)
            prob.model.add_objective('topo.ks_nodal_stress', scaler=1.0)
            prob.model.add_constraint('topo.m', upper=full_mass*in_design_mass)

            # Set up problem info for output
            prob_info = 'mass-{:.2f}'.format(in_design_mass)

        elif in_opt_problem == 'mass_min_stress_constr':

            # Check inputs
            if in_design_stress is None:
                raise ValueError("\ndesign_stress must be specified!\n")

            if in_stress_as_fraction:

                analysis = PlaneStressAnalysis(conn, dof, X, force, r0,
                    C, density, in_qval, in_epsilon, in_ks_parameter)
                x_full = np.ones(nnodes)
                stress_full = np.max(analysis.nodal_stress(x_full))
                design_stress = in_design_stress*stress_full

            else:
                design_stress = in_design_stress

            # Save to pkl file for postprocess' purpose
            prob_pkl['design_stress'] = design_stress

            # Create analysis object
            analysis = PlaneStressAnalysis(conn, dof, X, force, r0,
                    C, density, in_qval, in_epsilon, in_ks_parameter,
                    compute_stress=True, compute_mass=True, design_stress=design_stress)

            # Define openmdao problem
            prob.model.add_subsystem('topo', analysis)
            prob.model.connect('indeps.x', 'topo.x')
            prob.model.add_design_var('indeps.x', lower=1e-3, upper=1.0)
            prob.model.add_objective('topo.m', scaler=1.0)
            prob.model.add_constraint('topo.ks_nodal_stress', upper=1.0)

            # Set up problem info for output
            prob_info = 'stress-{:.2e}'.format(design_stress)

        # Set up output name
        outputname = '{:s}-{:s}-{:s}-{:s}-q{:.1f}-e{:.1f}-ks{:.1f}-iter{:d}'.format(
            in_opt_problem, optimizer_info, prob_info, prob_name,
            in_qval, in_epsilon, in_ks_parameter, in_max_iter)
        if in_outdir is not None:
            outputname = in_outdir +'/' + outputname

        # Setup optimizer
        if in_optimizer == 'ParOpt':
            if in_ParOpt_use_filter:
                tr_accept_step_strategy = 'filter_method'
            else:
                tr_accept_step_strategy = 'penalty_method'
            prob.driver = ParOptDriver()
            options = {
                'algorithm': 'tr',
                'output_level':0,
                'norm_type': 'l1',
                'tr_init_size': 0.05,
                'tr_min_size': 1e-3,
                'tr_max_size': 10.0,
                'tr_eta': 0.25,
                'tr_infeas_tol': 1e-6,
                'tr_l1_tol': 0.0,
                'tr_linfty_tol': 0.0,
                'tr_adaptive_gamma_update':
                    in_ParOpt_use_adaptive_gamma_update,
                'tr_accept_step_strategy':
                    tr_accept_step_strategy,
                'filter_sufficient_reduction': True,
                'filter_has_feas_restore_phase': True,
                'tr_use_soc': in_ParOpt_use_soc,
                'tr_max_iterations': in_max_iter,
                'output_file': outputname+'.out',
                'tr_output_file': outputname+'.tr',
                'penalty_gamma': ParOpt_penalty_gamma,
                'qn_subspace_size': qn_subspace_size,
                'qn_type': 'bfgs',
                'qn_diag_type': ParOpt_qn_diag_type,
                'abs_res_tol': 1e-8,
                'starting_point_strategy': 'affine_step',
                'barrier_strategy': 'mehrotra_predictor_corrector',
                'tr_steering_barrier_strategy':
                    'mehrotra_predictor_corrector',
                'tr_steering_starting_point_strategy': 'affine_step',
                'use_line_search': False}
            for key in options:
                prob.driver.options[key] = options[key]

        elif in_optimizer == 'IPOPT':
            prob.driver = om.pyOptSparseDriver()
            prob.driver.options['optimizer'] = 'IPOPT'
            prob.driver.opt_settings['max_iter'] = in_max_iter
            prob.driver.opt_settings['tol'] = 1e-10
            prob.driver.opt_settings['constr_viol_tol'] = 1e-10
            prob.driver.opt_settings['dual_inf_tol'] = 1e-10
            prob.driver.opt_settings['output_file'] = outputname+'.out'

        elif in_optimizer == 'SNOPT':
            prob.driver = om.pyOptSparseDriver()
            prob.driver.options['optimizer'] = 'SNOPT'
            prob.driver.opt_settings['Iterations limit'] = 9999999999999
            prob.driver.opt_settings['Major feasibility tolerance'] = 1e-10
            prob.driver.opt_settings['Major optimality tolerance'] = 1e-10
            prob.driver.opt_settings['Major iterations limit'] = in_max_iter
            prob.driver.opt_settings['Summary file'] = outputname+'.out'
            prob.driver.opt_settings['Print file'] = outputname+'_print.out'
            prob.driver.opt_settings['Major print level'] = 1
            prob.driver.opt_settings['Minor print level'] = 0

        # Run optimization and time it
        prob.setup()

        t_start = timeit.default_timer()
        prob.run_driver()
        t_end = timeit.default_timer()

        # Get design
        x_opt = prob.get_val('indeps.x')

        # Update starting point for next continuation step
        x_temp = x_opt

        # Plot history
        if in_optimizer == 'ParOpt':
            try:
                from paropt import plot_history
            except:
                print("\n[Warning] cannot import plot_history from paropt,",
                    "no convergence history plot generated!\n")
            else:
                tr_name = outputname + '.tr'
                plot_history(tr_name, savefig=True)

        if in_optimizer == 'IPOPT':
            try:
                from plane_stress.utils import ipopt_plot
            except:
                print("\n[Warning] cannot import ipopt_plot from plane_stress.utils,",
                    "no convergence history plot generated!\n")
            else:
                out_name = outputname + '.out'
                ipopt_plot(out_name, savefig=True)

        if in_optimizer == 'SNOPT':
            try:
                from plane_stress.utils import snopt_plot
            except:
                print("\n[Warning] cannot import snopt_plot from plane_stress.utils,",
                    "no convergence history plot generated!\n")
            else:
                out_name = outputname + '.out'
                snopt_plot(out_name, savefig=True)

        index += 1


    # Plot design and stress contour for final design
    topo_name = outputname + '_topo.png'
    stress_name = outputname + '_stress.png'
    analysis.plot_topology(x_opt, savefig=True, filename=topo_name)
    analysis.plot_stress(x_opt, savefig=True, filename=stress_name)

    # Save extra information to pkl file
    prob_pkl['x'] = x_opt
    prob_pkl['opt_time'] = '{:.2e} s'.format(t_end - t_start)
    prob_pkl['qval'] = in_qval
    prob_pkl['epsilon'] = in_epsilon
    prob_pkl['ks_parameter'] = in_ks_parameter

    # Output pkl file
    picklename = outputname+'.pkl'
    with open(picklename, 'wb') as pklfile:
        pickle.dump(prob_pkl, pklfile)


if __name__ == '__main__':

    # Set up parser
    p = argparse.ArgumentParser()
    p.add_argument('picklename', type=str,
        help='pickle file that defines the analysis problem')
    p.add_argument('opt_problem', type=str, choices=[
        'comp_min_mass_constr', 'comp_min_massfreq_constr',
        'comp_min_massstress_constr', 'comp_min_massfreqstress_constr',
        'stress_min_mass_constr', 'mass_min_stress_constr'])
    p.add_argument('optimizer', type=str, choices=['ParOpt', 'SNOPT', 'IPOPT'])
    p.add_argument('--qval', nargs='*', type=float, default=[5.0])
    p.add_argument('--epsilon', type=float, default=0.1)
    p.add_argument('--ks_parameter', type=float, default=100.0)
    p.add_argument('--design_mass', type=float, default=None,
        help='normalized against full material design')
    p.add_argument('--design_freq', type=float, default=None)
    p.add_argument('--design_stress', type=float, default=None)
    p.add_argument('--stress_as_fraction', action='store_true')
    p.add_argument('--max_iter', nargs='*', type=int, default=[2])
    p.add_argument('--ParOpt_use_adaptive_gamma_update', action='store_true')
    p.add_argument('--ParOpt_use_filter', action='store_true')
    p.add_argument('--ParOpt_use_soc', action='store_true')
    p.add_argument('--info', type=str, default=None)
    p.add_argument('--outdir', type=str, default=None)
    args = p.parse_args()

    optimize(args.picklename, args.opt_problem, args.optimizer,
             args.qval, args.epsilon, args.ks_parameter, args.design_mass,
             args.design_freq, args.design_stress, args.stress_as_fraction,
             args.max_iter, args.ParOpt_use_adaptive_gamma_update,
             args.ParOpt_use_filter, args.ParOpt_use_soc, args.info, args.outdir)
