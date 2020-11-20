#!/usr/bin/env python3

"""
This script generates performance profiler(s) given a batch of cases.
"""

import numpy as np
import argparse
import os
import pickle
import matplotlib.pyplot as plt

# Define colors
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
        '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
        '#bcbd22', '#17becf']

# Name of all optimizers
optimizers = ['ParOpt',
              'ParOptFilter',
              'ParOptFilterSoc',
              'ParOptAdapt',
              'IPOPT',
              'SNOPT',
              ]

# Parser
p = argparse.ArgumentParser()
p.add_argument('--result_folder', type=str, default='results')
p.add_argument('--opt_problem', nargs='*', type=str, default=None, choices=[
    'all', 'comp_min_mass_constr', 'comp_min_massfreq_constr',
    'comp_min_massstress_constr', 'comp_min_massfreqstress_constr',
    'stress_min_mass_constr', 'mass_min_stress_constr'])
p.add_argument('--metric', type=str, default='obj', choices=['obj', 'time'])
p.add_argument('--infeas_tol', type=float, default=1e-4)
p.add_argument('--showfig', action='store_true')
args = p.parse_args()

# Get problem list
if 'all' in args.opt_problem:
    opt_problems = ['comp_min_mass_constr', 'comp_min_massfreq_constr',
        'comp_min_massstress_constr' , 'mass_min_stress_constr']
else:
    opt_problems = args.opt_problem

# Input value check
if args.opt_problem is None:
    raise ValueError("\n At least one problem need to be specified using --opt_problem!")

# Get which performance metric to use
if args.metric == 'obj':
    metric = 'obj'
elif args.metric == 'time':
    metric = 'opt_time'

# Store all folder names
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

# Get total number of different optimization problems
nprobs = len(case_folders)

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
for case_folder in case_folders:

    # Loop over each files
    for f in os.listdir(args.result_folder+'/'+case_folder):

        # We only need result pkl files
        if os.path.splitext(f)[1] == '.pkl':

            # Get optimizer related to this pkl
            optimizer_used = f.split('-')[1]

            with open('{:s}/{:s}/{:s}'.format(
                args.result_folder, case_folder, f), 'rb') as pklfile:
                prob_pkl = pickle.load(pklfile)

                # Optimization time is stored as string in pkl, we convert it to time
                # in case time is used as metric
                prob_pkl['opt_time'] = float(prob_pkl['opt_time'][:-1])

                # Load obj/time
                try:
                    metricval = float(prob_pkl[metric])
                    optimizer_data[optimizer_used]['metric'][prob_index] = metricval
                except:
                    pass

                # Load maximum normalized constraint violation
                # TODO: fix constraint in optimize
                try:
                    # We read constraints from text output file
                    if 'ParOpt' in optimizer:
                        textout = '{:s}/{:s}/{:s}.tr'.format(args.result_folder, case_folder,
                            os.path.splitext(f)[0])
                        with open(textout, 'r') as ff:
                            for line in ff:
                                pass
                            infeas = float(line.split()[2])
                    else:
                        textout = '{:s}/{:s}/{:s}.out'.format(args.result_folder, case_folder,
                            os.path.splitext(f)[0])
                        if optimizer == 'IPOPT':
                            with open(textout, 'r') as ff:
                                for line in ff:
                                    if 'Constraint violation....:' in line:
                                        infeas = float(line.split()[-1])
                                        break
                        elif optimizer == 'SNOPT':
                            with open(textout, 'r') as ff:
                                for line in ff:
                                    if 'Nonlinear constraint violn' in line:
                                        infeas = float(line.split()[-1])
                                        break
                    optimizer_data[optimizer_used]['infeas'][prob_index] = infeas
                    if optimizer_data[optimizer_used]['infeas'][prob_index] > args.infeas_tol:
                        optimizer_data[optimizer_used]['metric'][prob_index] = None
                except:
                    optimizer_data[optimizer_used]['metric'][prob_index] = None

                # print metric, infeas and file name to output for checking purpose
                try:
                    print("{:s}/{:s}/{:s}\nmetric: {:e}  infeas: {:e}".format(
                        args.result_folder, case_folder, f, metricval, infeas))
                except:
                    pass

    # Now normalize the metric against the best one
    metrics = [optimizer_data[key]['metric'][prob_index] for key in optimizers]
    values = [val for val in metrics if type(val) is float]

    for optimizer in optimizers:
        if type(optimizer_data[optimizer]['metric'][prob_index]) is float:
            best = min(values)
            optimizer_data[optimizer]['metric'][prob_index] /= best

    prob_index += 1

# Now we can plot the profiler
fig, ax1 = plt.subplots()

# One curve for each optimizer
x_lim = 0.0
for optimizer in optimizers:
    sorted_metrics = sorted([optimizer_data[optimizer]['metric'][i] for i in range(nprobs)
        if type(optimizer_data[optimizer]['metric'][i]) is float])
    percentiles = np.linspace(1, len(sorted_metrics), len(sorted_metrics))
    percentiles /= nprobs

    optimizer_data[optimizer]['sorted_metrics'] = sorted_metrics
    optimizer_data[optimizer]['percentiles'] = percentiles

    try:
        if max(sorted_metrics) > x_lim:
            x_lim = max(sorted_metrics)
    except:
        pass

index = 0
for optimizer in optimizers:
    optimizer_data[optimizer]['sorted_metrics'].append(x_lim)
    optimizer_data[optimizer]['percentiles'] = np.append(optimizer_data[optimizer]['percentiles'],
        optimizer_data[optimizer]['percentiles'][-1])

    ax1.step(optimizer_data[optimizer]['sorted_metrics'],
             optimizer_data[optimizer]['percentiles'],
             label=optimizer, color=colors[index])

    index += 1

ax1.set_xlabel('Normalized objective')
ax1.set_ylabel('Percentage of cases')
plt.title('Performance profiler')
plt.legend()

if args.showfig:
    plt.show()
else:
    name = 'profiler'
    for prob in opt_problems:
        name += '-'+prob
    name += '.png'
    plt.savefig(name)
    plt.close()