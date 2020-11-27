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
              'IPOPT',
              'SNOPT',
              ]

legends = {
    'ParOpt': r'ParOpt S$\ell_{1}$QP',
    'ParOptAdapt': r'ParOpt S$\ell_{1}$QP w/ adaptive',
    'ParOptFilter': r'ParOpt filterSQP',
    'ParOptFilterSoc': r'ParOpt filterSQP w/ SOC',
    'IPOPT': r'IPOPT',
    'SNOPT': r'SNOPT',

}

# Parser
p = argparse.ArgumentParser()
p.add_argument('--result_folder', type=str, default='results')
p.add_argument('--opt_problem', nargs='*', type=str, default=['all'], choices=[
    'all', 'comp_min_mass_constr', 'comp_min_massfreq_constr',
    'comp_min_massstress_constr', 'comp_min_massfreqstress_constr',
    'stress_min_mass_constr', 'mass_min_stress_constr'])
p.add_argument('--metric', type=str, default='obj', choices=['obj', 'time', 'discreteness'])
p.add_argument('--infeas_tol', type=float, default=1e-4)
p.add_argument('--metric_limit', type=float, default=2.0,
    help='cases with metrics larger than this limit are considered failed')
p.add_argument('--plot', action='store_true')
args = p.parse_args()

# Get problem list
if 'all' in args.opt_problem:
    opt_problems = ['comp_min_mass_constr', 'comp_min_massfreq_constr',
        'comp_min_massstress_constr' , 'mass_min_stress_constr']
else:
    opt_problems = args.opt_problem

# Get which performance metric to use
if args.metric == 'obj':
    metric = 'obj'
elif args.metric == 'time':
    metric = 'opt_time'

metric_xmin = 0.95
metric_xmax = args.metric_limit

if args.metric == 'discreteness':
    metric_xmin = -0.01
    metric_xmax = 0.25

# Store all folder names in results/ folder
case_folders = []

# Store all domain names, which is the same as the name of mesh pickle
case_domains = []

# Loop over all folders in results/
for opt_problem in opt_problems:
    for folder_in_results in os.listdir(args.result_folder):
        if opt_problem in folder_in_results:
            case_folders.append(folder_in_results)
            if folder_in_results.replace(opt_problem+'-','') not in case_domains:
                case_domains.append(folder_in_results.replace(opt_problem+'-',''))

# Since each folder contains all cases for one specific optimization problem,
# we get total number of different optimization problems by the following
nprobs = len(case_folders)

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
for case_folder in case_folders:

    # Store all pkl files
    pkl_files = []

    # Loop over files in each case folder
    for case_file in os.listdir(args.result_folder+'/'+case_folder):

        # Gget name of all pkl files:
        if os.path.splitext(case_file)[1] == '.pkl':
            pkl_files.append(case_file)

    # If number of pkls is less than number of optimizers, we skip this optimization
    # problem
    if len(pkl_files) < len(optimizers):
        print('[Warning] Missing cases in the following folder!\n{:s}/{:s}'.format(
            args.result_folder, case_folder))
        continue
    else:
        nvalids += 1

    # Next, loop over each pkl to get metric and infeas
    for pkl_file in pkl_files:

        # Get optimizer name
        optimizer = pkl_file.split('-')[1]

        # Load pkl
        with open('{:s}/{:s}/{:s}'.format(
            args.result_folder, case_folder, pkl_file), 'rb') as p:
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
                        args.result_folder, case_folder, pkl_file))
                    print("Error message:", e)
            else:
                try:
                    optimizer_data[optimizer]['metric'][prob_index] = float(pkl_dict[metric])

                except Exception as e:
                    print("\n[Warning] cannot load pkl file:\n{:s}/{:s}/{:s}".format(
                        args.result_folder, case_folder, pkl_file))
                    print("Error message:", e)

            # Load maximum normalized constraint violation
            # We read constraints from text output file
            if 'ParOpt' in optimizer:
                textout = '{:s}/{:s}/{:s}.tr'.format(args.result_folder, case_folder,
                    os.path.splitext(pkl_file)[0])
                with open(textout, 'r') as ff:
                    for line in ff:
                        pass
                    infeas = float(line.split()[2])
            else:
                textout = '{:s}/{:s}/{:s}.out'.format(args.result_folder, case_folder,
                    os.path.splitext(pkl_file)[0])
                if optimizer == 'IPOPT':
                    with open(textout, 'r') as ff:
                        for line in ff:
                            if 'Constraint violation' in line:
                                infeas = float(line.split()[-1])
                                break
                elif optimizer == 'SNOPT':
                    with open(textout, 'r') as ff:
                        for line in ff:
                            if 'Nonlinear constraint violn' in line:
                                infeas = float(line.split()[-1])
                                break

            optimizer_data[optimizer]['infeas'][prob_index] = infeas

            print('-----------------------------------------------------------')
            print('[pkl] {:s}/{:s}/{:s}'.format(args.result_folder, case_folder, pkl_file))
            print('[obj] {:f}'.format(optimizer_data[optimizer]['metric'][prob_index]))
            print('[con] {:f}'.format(optimizer_data[optimizer]['infeas'][prob_index]))

            if optimizer_data[optimizer]['infeas'][prob_index] > args.infeas_tol:
                optimizer_data[optimizer]['metric'][prob_index] = None
                ninfeas += 1

    # Now we've get all data needed in this problem folder
    # Next, normalize the metric against the best one
    metrics = [ optimizer_data[key]['metric'][prob_index] for key in optimizers if
        optimizer_data[key]['metric'][prob_index] is not None ]

    if metrics:
        # compute best over all optimziers if metrics is not empty
        best = min(metrics)
    else:
        # Otherwise, give a dummy quantity to best
        best = 1.0
        print('[Warning] No feasible cases in the following folder!\n{:s}/{:s}'.format(
            args.result_folder, case_folder))
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
print('number of individual infeasible cases                  :  {:d}'.format(ninfeas))
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
    ax.step(sorted_metrics, percentiles, label=legends[optimizer])

    index += 1

if len(opt_problems) == 4:
    title = 'Performance Profiler on Full 2D Problem Set'
    name = 'profiler-{:s}-all-tol-{:.0e}-ninfeas-{:d}.pdf'.format(args.metric, args.infeas_tol, ninfeas)

elif len(opt_problems) == 1:
    title = 'Performance Profiler on Problem Set: {:s}'.format(opt_problems[0])
    name = 'profiler-{:s}-{:s}-tol-{:.0e}-ninfeas-{:d}.pdf'.format(args.metric, opt_problems[0], args.infeas_tol, ninfeas)
else:
    title = 'Performance Profiler'
    name = 'profiler'
    for prob in opt_problems:
        name += '-'+prob
    name += '-tol-{:.0e}-ninfeas-{:d}.pdf'.format(args.infeas_tol, ninfeas)

ax.set_xlim([metric_xmin, metric_xmax])
ax.set_ylim([0.0, 1.01])

if args.metric == 'obj':
    xlabel = 'Normalized objective'
elif args.metric == 'time':
    xlabel = r'Normalized optimization time \tau'
elif args.metric == 'discreteness':
    xlabel = 'Average discreteness'
ax.set_xlabel(xlabel)
ax.set_ylabel('Percentage of cases')
plt.legend()

if args.plot:
    plt.show()
else:
    plt.savefig(name)
    plt.close()