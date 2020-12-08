#!/usr/bin/env python3

"""
This script generates performance profiler(s) given a batch of cases.

First, we define some terminologies. Usually, a batch of cases consist
of multiple dimensions, e.g. meshsize, meshtype, domain shape, etc. Also
for a same physical problem with identical geometry, mesh and boundary
condition, we could formulate different optimization problems with different
objective and constraint(s) of interest.

Thus, a case can vary in the following three categories:
    - physical problem
    - optimization problem
    - optimizer

where,
    1. physical problem defines meshtype (structured, unstructured),
       problem type (cantilever, michell, MBB, lbracket), Aspect
       Ratio (1.0, 2.0, 3.0) and whether there is a hole in domain or not.
       In our implementation, we use preprocessor script to generate physical
       problem, and save information in python pickle file *.pkl.

    2. optimization problem contains physical problem too, besides, it also
       defines objective and constraint(s), we
       currently implemented 6 optimization problems:
       - comp_min_mass_constr
       - comp_min_massfreq_constr
       - comp_min_massstress_constr
       - comp_min_massfreqstress_constr
       - stress_min_mass_constr
       - mass_min_stress_constr

    3. optimizer decides optimization settings, which algorithm to use, etc.
"""

import numpy as np
import argparse
import os
import pickle
import matplotlib.pyplot as plt
from pprint import pprint

# Name of all optimizers
optimizers = ['ParOpt',
              'ParOptAdapt',
              'ParOptFilter',
              'ParOptFilterSoc',
              'ParOptQn',
              'ParOptAdaptQn',
              'ParOptFilterQn',
              'ParOptFilterSocQn',
              'IPOPT',
              'IPOPTtuned',
              'SNOPT',
              'SNOPTtuned',
            #   'SNOPTtunedscaled'
              ]

optimizers_no_qn = ['ParOpt',
              'ParOptAdapt',
              'ParOptFilter',
              'ParOptFilterSoc',
              'IPOPT',
              'IPOPTtuned',
              'SNOPT',
              'SNOPTtuned',
            #   'SNOPTtunedscaled'
              ]

legends = {
    'ParOpt': r'ParOpt S$\ell_{1}$QP',
    'ParOptAdapt': r'ParOpt S$\ell_{1}$QP w/ adaptive',
    'ParOptFilter': r'ParOpt filterSQP',
    'ParOptFilterSoc': r'ParOpt filterSQP w/ SOC',
    'ParOptQn': r'ParOpt S$\ell_{1}$QP w/ qn',
    'ParOptAdaptQn': r'ParOpt S$\ell_{1}$QP w/ adaptive and qn',
    'ParOptFilterQn': r'ParOpt filterSQP w/ qn',
    'ParOptFilterSocQn': r'ParOpt filterSQP w/ SOC and qn',
    'IPOPT': r'IPOPT',
    'SNOPT': r'SNOPT',
    'IPOPTtuned': r'IPOPT tuned',
    'SNOPTtuned': r'SNOPT tuned',
    'SNOPTtunedscaled': r'SNOPT scaled'
}

colors = {
    'ParOpt': '#DE425B',
    'ParOptAdapt': '#FF6361',
    'ParOptFilter': '#00876C',
    'ParOptFilterSoc': '#89BF77',
    'ParOptQn': '#FFCAB6',
    'ParOptAdaptQn': '#FFA600',
    'ParOptFilterQn': '#58508D',
    'ParOptFilterSocQn': '#BC5090',
    'IPOPT': '#7A4EFE',
    'IPOPTtuned': '#00ADff',
    'SNOPT': '#2e2e2e',
    'SNOPTtuned': '#898989',
    'SNOPTtunedscaled': 'black'
}

# Parser
p = argparse.ArgumentParser()
p.add_argument('--result_folder', nargs='*', type=str, default=['results'])
p.add_argument('--opt_problem', nargs='*', type=str, default=['all'], choices=[
    'all', 'allcomp', 'comp_min_mass_constr', 'comp_min_massfreq_constr',
    'comp_min_massstress_constr', 'comp_min_massfreqstress_constr',
    'stress_min_mass_constr', 'mass_min_stress_constr'])
p.add_argument('--metric', type=str, default='obj', choices=['obj', 'time', 'discreteness'])
p.add_argument('--infeas_tol', type=float, default=1e-4)
p.add_argument('--metric_limit', type=float, default=2.0,
    help='cases with metrics larger than this limit are considered failed')
p.add_argument('--separate_legend', action='store_true')
p.add_argument('--plot', action='store_true')
args = p.parse_args()

# Get problem list
if 'all' in args.opt_problem:
    opt_problems = ['comp_min_mass_constr', 'comp_min_massfreq_constr',
        'comp_min_massstress_constr' , 'mass_min_stress_constr']
elif 'allcomp' in args.opt_problem:
     opt_problems = ['comp_min_mass_constr', 'comp_min_massfreq_constr',
        'comp_min_massstress_constr']
else:
    opt_problems = args.opt_problem

# we don't plot profiles with qn correction if there is non-compliance
# cases in the problem set
if ('stress_min_mass_constr' in opt_problems or
    'mass_min_stress_constr' in opt_problems):
    optimizers = optimizers_no_qn

# Get which performance metric to use
if args.metric == 'obj':
    metric = 'obj'
elif args.metric == 'time':
    metric = 'opt_time'

metric_xmin = 1.0 - 0.05*(args.metric_limit - 1.0)
metric_xmax = args.metric_limit

if args.metric == 'discreteness':
    metric_xmin = -0.01
    metric_xmax = 0.25

# Store all folder names in results/ folder
case_folders = dict()

# Store all domain names, which is the same as the name of mesh pickle
case_domains = []
nprobs = 0
# Loop over all folders in results/
for result_folder in args.result_folder:
    case_folders[result_folder] = []
    for opt_problem in opt_problems:
        for folder_in_results in os.listdir(result_folder):
            if opt_problem in folder_in_results:
                case_folders[result_folder].append(folder_in_results)
                if folder_in_results.replace(opt_problem+'-','') not in case_domains:
                    case_domains.append(folder_in_results.replace(opt_problem+'-',''))

    # Since each folder contains all cases for one specific optimization problem,
    # we get total number of different optimization problems by the following
    nprobs += len(case_folders[result_folder])

# We also want to keep track of number of valid problems, nvalids <= nprobs
nvalids = 0

# Use dictionary to store objective and max infeasibility for each optimizer
optimizer_data = dict()

# Initialize dictionary
for optimizer in optimizers:
    optimizer_data[optimizer] = dict()

    # Initialize obj and infeas lists, note that
    # None indicating no results found yet
    optimizer_data[optimizer]['metric'] = [None] * nprobs
    optimizer_data[optimizer]['infeas'] = [None] * nprobs

# Get objectives and max constraint violation for each problem for each optimizer
prob_index = 0
ninfeas = 0
ninfeas_all_optimizer = 0
nplotted = 0
n_obj_near_zero = 0
for result_folder in args.result_folder:
    for case_folder in case_folders[result_folder]:

        # Store all pkl files
        pkl_files = []

        # Loop over files in each case folder
        for case_file in os.listdir(result_folder+'/'+case_folder):

            # Gget name of all pkl files:
            if os.path.splitext(case_file)[1] == '.pkl':

                # Get optimizer name
                _optimizer = case_file.split('-')[1]
                if _optimizer in optimizers:
                    pkl_files.append(case_file)

        # If number of pkls is not equal to number of optimizers, we skip this optimization
        # problem
        if len(pkl_files) != len(optimizers):
            print('[Warning] Missing cases in the following folder!\n{:s}/{:s}'.format(
                result_folder, case_folder))
            continue
        else:
            nvalids += 1

        # Next, loop over each pkl to get metric and infeas
        for pkl_file in pkl_files:

            # Get optimizer name
            optimizer = pkl_file.split('-')[1]

            # We skip the pickles using the optimizer that is not in the optimizer list
            if optimizer not in optimizers:
                continue

            # Load pkl
            with open('{:s}/{:s}/{:s}'.format(
                result_folder, case_folder, pkl_file), 'rb') as p:
                pkl_dict = pickle.load(p)

                # Optimization time is stored as string in pkl, we convert it to time
                # in case time is used as metric
                pkl_dict['opt_time'] = float(pkl_dict['opt_time'][:-1])

                # Get metric
                if args.metric == 'discreteness':
                    try:
                        x = pkl_dict['x']
                        optimizer_data[optimizer]['metric'][prob_index] = np.dot(x, 1-x)/len(x)
                    except Exception as e:
                        print("\n[Warning] cannot load pkl file:\n{:s}/{:s}/{:s}".format(
                            result_folder, case_folder, pkl_file))
                        print("Error message:", e)
                else:
                    try:
                        optimizer_data[optimizer]['metric'][prob_index] = float(pkl_dict[metric])

                    except Exception as e:
                        print("\n[Warning] cannot load pkl file:\n{:s}/{:s}/{:s}".format(
                            result_folder, case_folder, pkl_file))
                        print("Error message:", e)

                # Load maximum normalized constraint violation
                # We read constraints from text output file
                if 'ParOpt' in optimizer:
                    textout = '{:s}/{:s}/{:s}.tr'.format(result_folder, case_folder,
                        os.path.splitext(pkl_file)[0])
                    with open(textout, 'r') as ff:
                        for line in ff:
                            pass
                        infeas = float(line.split()[2])
                else:
                    textout = '{:s}/{:s}/{:s}.out'.format(result_folder, case_folder,
                        os.path.splitext(pkl_file)[0])
                    if 'IPOPT' in optimizer:
                        with open(textout, 'r') as ff:
                            for line in ff:
                                if 'Constraint violation' in line:
                                    infeas = float(line.split()[-1])
                                    break
                    elif 'SNOPT' in optimizer:
                        with open(textout, 'r') as ff:
                            for line in ff:
                                if 'Nonlinear constraint violn' in line:
                                    infeas = float(line.split()[-1])
                                    break

                optimizer_data[optimizer]['infeas'][prob_index] = infeas

                # print('-----------------------------------------------------------')
                # print('[pkl] {:s}/{:s}/{:s}'.format(result_folder, case_folder, pkl_file))
                # print('[obj] {:f}'.format(optimizer_data[optimizer]['metric'][prob_index]))
                # print('[con] {:f}'.format(optimizer_data[optimizer]['infeas'][prob_index]))

                if optimizer_data[optimizer]['infeas'][prob_index] > args.infeas_tol:
                    optimizer_data[optimizer]['metric'][prob_index] = None
                    ninfeas += 1
                else:
                    nplotted += 1

        # Now we've get all data needed in this problem folder
        # Next, normalize the metric against the best one
        metrics = [ optimizer_data[key]['metric'][prob_index] for key in optimizers if
            optimizer_data[key]['metric'][prob_index] is not None ]

        if metrics:
            # compute best over all optimziers if metrics is not empty
            best = min(metrics)
            if best < 1e-15:
                best = 1e-15
                print('[Warning] objective near zero found!')
                n_obj_near_zero += 1
        else:
            # Otherwise, give a dummy quantity to best
            best = 1.0
            print('[Warning] No feasible cases in the following folder!\n{:s}/{:s}'.format(
                result_folder, case_folder))
            ninfeas_all_optimizer += 1

        # We don't normalize discreteness against the best because it is already dimentionless
        if args.metric == 'time' or args.metric == 'obj':
            for optimizer in optimizers:
                if optimizer_data[optimizer]['metric'][prob_index] is not None:
                    optimizer_data[optimizer]['metric'][prob_index] /= best

        prob_index += 1

# Print out summary
print('\n')
print('                -------SUMMARY-------')
print('number of optimizers                                   :  {:d}'.format(len(optimizers)))
print('number of specified optimization problems              :  {:d}'.format(nprobs))
print('number of valid problems (every optimizer has a result):  {:d}'.format(nvalids))
print('number of problems that have no feasible cases at all  :  {:d}'.format(ninfeas_all_optimizer))
print('number of specified optimization cases                 :  {:d}'.format(len(optimizers)*nprobs))
print('number of plotted optimization cases                   :  {:d}/{:d} ({:.2f}%)'.format(nplotted, len(optimizers)*nprobs, 100*nplotted/len(optimizers)/nprobs))
print('number of individual infeasible cases                  :  {:d}/{:d} ({:.2f}%)'.format(ninfeas, len(optimizers)*nprobs, 100*ninfeas/len(optimizers)/nprobs))
print('number of best objectives near zero                    :  {:d}'.format(n_obj_near_zero))
print('\n')

# Set up plotting environment
mpl_style_path = os.path.dirname(os.path.realpath(__file__)) + '/paper.mplstyle'
plt.style.use(mpl_style_path)
fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8))

index = 0
for optimizer in optimizers:

    # We sort all valid normalized metrics for current optimizer
    sorted_metrics = sorted([ optimizer_data[optimizer]['metric'][i] for i in range(nprobs)
        if optimizer_data[optimizer]['metric'][i] is not None ])

    # We compute percentile based on actual total number of cases
    percentiles = [ (sorted_metrics.index(i) + 1 ) / nvalids for i in sorted_metrics ]

    # Append one more entry so that we will have a flat curve till right end
    if not sorted_metrics or sorted_metrics[-1] < metric_xmax:
        sorted_metrics.append(metric_xmax)
        if percentiles:
            percentiles.append(percentiles[-1])
        else:
            percentiles.append(0.0)

    # Insert (1.0, 0.0) to beginning of lists so that we will have all curves
    # starting from (1.0, 0.0)
    # sorted_metrics.insert(0, 1.0)
    # percentiles.insert(0, 0.0)

    # Plot
    ax.step(sorted_metrics, percentiles, label=legends[optimizer], color=colors[optimizer])

    index += 1

# if len(opt_problems) == 4:
#     title = 'Performance Profiler on Full Problem Set'
#     name = 'profiler-{:s}-all-tol-{:.0e}-ninfeas-{:d}.pdf'.format(args.metric, args.infeas_tol, ninfeas)

# if len(opt_problems) == 1:
#     title = 'Performance Profiler on Problem Set: {:s}'.format(opt_problems[0])
#     name = 'profiler-{:s}-{:s}-tol-{:.0e}-ninfeas-{:d}.pdf'.format(args.metric, opt_problems[0], args.infeas_tol, ninfeas)
# else:
title = 'Performance Profiler'
name = 'profiler-{:s}'.format(args.metric)
for prob in opt_problems:
    name += '-'+prob
name += '-tol-{:.0e}-ninfeas-{:d}-nplotted-{:d}-nfolders-{:d}.pdf'.format(args.infeas_tol, ninfeas, nplotted, len(args.result_folder))

ax.set_xlim([metric_xmin, metric_xmax])
ax.set_ylim([0.0, 1.01])

if args.metric == 'obj':
    xlabel = 'Normalized objective'
elif args.metric == 'time':
    xlabel = r'Normalized optimization time \tau'
elif args.metric == 'discreteness':
    xlabel = 'Average discreteness'
ax.set_xlabel(xlabel)
ax.set_ylabel('Fraction of cases')

if not args.separate_legend:
    plt.legend()

if args.plot:
    plt.show()
else:
    plt.savefig(name)
    legendfig = plt.figure()

    # Save legend separately
    if args.separate_legend:
        plt.figlegend(*ax.get_legend_handles_labels(), loc='upper left', ncol=3)
        plt.savefig('legend-{:d}-optimizers.pdf'.format(len(optimizers)))

    # close the plt
    plt.close()