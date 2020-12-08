#!/usr/bin/env python3

import os
import argparse

p = argparse.ArgumentParser()
p.add_argument('--detail', action='store_true')
args = p.parse_args()

optimizer_set_basic = ['ParOpt', 'ParOptAdapt', 'ParOptFilter',
    'ParOptFilterSoc', 'SNOPT', 'IPOPT']
optimizer_set_qn = ['ParOptQn', 'ParOptAdaptQn', 'ParOptFilterQn',
    'ParOptFilterSocQn']
optimizer_set_tune = ['IPOPTtuned', 'SNOPTtuned']
optimizer_set_scale = ['SNOPTtunedscaled']

# Main result folder
main_folder = 'results'

# Enumate all optimization problems:
opt_problem_counters = {
    'comp_min_mass_constr': {
        'n_opt_basic': 0,
        'n_opt_basic_target': 0,
        'n_opt_qn': 0,
        'n_opt_qn_target': 0,
        'n_opt_tune': 0,
        'n_opt_tune_target': 0,
        'n_opt_scale': 0,
        'n_opt_scale_target': 0,
        'N_OPT_BASIC' : len(optimizer_set_basic),
        'N_OPT_QN': len(optimizer_set_qn),
        'N_OPT_TUNE': len(optimizer_set_tune),
        'N_OPT_SCALE': len(optimizer_set_scale),
        'N_OPT': len(optimizer_set_basic) + len(optimizer_set_qn) + len(optimizer_set_tune) + len(optimizer_set_scale)
    },
    'comp_min_massfreq_constr': {
        'n_opt_basic': 0,
        'n_opt_basic_target': 0,
        'n_opt_qn': 0,
        'n_opt_qn_target': 0,
        'n_opt_tune': 0,
        'n_opt_tune_target': 0,
        'n_opt_scale': 0,
        'n_opt_scale_target': 0,
        'N_OPT_BASIC' : len(optimizer_set_basic),
        'N_OPT_QN': len(optimizer_set_qn),
        'N_OPT_TUNE': len(optimizer_set_tune),
        'N_OPT_SCALE': len(optimizer_set_scale),
        'N_OPT': len(optimizer_set_basic) + len(optimizer_set_qn) + len(optimizer_set_tune) + len(optimizer_set_scale)
    },
    'comp_min_massstress_constr': {
        'n_opt_basic': 0,
        'n_opt_basic_target': 0,
        'n_opt_qn': 0,
        'n_opt_qn_target': 0,
        'n_opt_tune': 0,
        'n_opt_tune_target': 0,
        'n_opt_scale': 0,
        'n_opt_scale_target': 0,
        'N_OPT_BASIC' : len(optimizer_set_basic),
        'N_OPT_QN': len(optimizer_set_qn),
        'N_OPT_TUNE': len(optimizer_set_tune),
        'N_OPT_SCALE': len(optimizer_set_scale),
        'N_OPT': len(optimizer_set_basic) + len(optimizer_set_qn) + len(optimizer_set_tune) + len(optimizer_set_scale)
    },
    'mass_min_stress_constr': {
        'n_opt_basic': 0,
        'n_opt_basic_target': 0,
        'n_opt_qn': 0,
        'n_opt_qn_target': 0,
        'n_opt_tune': 0,
        'n_opt_tune_target': 0,
        'n_opt_scale': 0,
        'n_opt_scale_target': 0,
        'N_OPT_BASIC' : len(optimizer_set_basic),
        'N_OPT_QN': None,
        'N_OPT_TUNE': len(optimizer_set_tune),
        'N_OPT_SCALE': len(optimizer_set_scale),
        'N_OPT': len(optimizer_set_basic) + len(optimizer_set_qn) + len(optimizer_set_tune) + len(optimizer_set_scale)
    }
}

# Initialize counter
counter_total = 0
counter_total_target = 0

# loop over all subfolders
subfolders = os.listdir(main_folder)
subfolders.sort()
for subfolder in subfolders:

    # Reset local counter
    counter_current = 0

    # Reset optimizer list
    local_optimizers = []

    # Get problem type
    opt_prob = subfolder.split('-')[0]

    # update target numbers
    if opt_problem_counters[opt_prob]['N_OPT_BASIC'] is not None:
        opt_problem_counters[opt_prob]['n_opt_basic_target'] += opt_problem_counters[opt_prob]['N_OPT_BASIC']
        counter_total_target += opt_problem_counters[opt_prob]['N_OPT_BASIC']

    if opt_problem_counters[opt_prob]['N_OPT_QN'] is not None:
        opt_problem_counters[opt_prob]['n_opt_qn_target'] += opt_problem_counters[opt_prob]['N_OPT_QN']
        counter_total_target += opt_problem_counters[opt_prob]['N_OPT_QN']

    if opt_problem_counters[opt_prob]['N_OPT_TUNE'] is not None:
        opt_problem_counters[opt_prob]['n_opt_tune_target'] += opt_problem_counters[opt_prob]['N_OPT_TUNE']
        counter_total_target += opt_problem_counters[opt_prob]['N_OPT_TUNE']

    if opt_problem_counters[opt_prob]['N_OPT_SCALE'] is not None:
        opt_problem_counters[opt_prob]['n_opt_scale_target'] += opt_problem_counters[opt_prob]['N_OPT_SCALE']
        counter_total_target += opt_problem_counters[opt_prob]['N_OPT_SCALE']

    # Loop over each pickle file in subfolders
    for f in os.listdir('{:s}/{:s}'.format(main_folder, subfolder)):
        if 'pkl' in f:

            opt_prob = f.split('-')[0]
            optimizer = f.split('-')[1]

            if optimizer in local_optimizers:
                print('\n[Warning] duplicate case found in {:s}/{:s} for optimizer {:s}'.format(
                    main_folder, subfolder, optimizer))

            else:
                local_optimizers.append(optimizer)
                if optimizer in optimizer_set_basic:
                    opt_problem_counters[opt_prob]['n_opt_basic'] += 1
                    counter_total += 1

                elif optimizer in optimizer_set_qn:
                    opt_problem_counters[opt_prob]['n_opt_qn'] += 1
                    counter_total += 1

                elif optimizer in optimizer_set_tune:
                    opt_problem_counters[opt_prob]['n_opt_tune'] += 1
                    counter_total += 1

                elif optimizer in optimizer_set_scale:
                    opt_problem_counters[opt_prob]['n_opt_scale'] += 1
                    counter_total += 1

                counter_current += 1

    # Print number of finished cases in each subfolder
    if args.detail:
        print('{:<90s}: {:>2d}/{:>2d}'.format(subfolder,
            counter_current, opt_problem_counters[opt_prob]['N_OPT']))

# Print summary
set_info = {
    'Basic set': {
        'n': 'n_opt_basic',
        'target': 'n_opt_basic_target'
    },
    'Qn set': {
        'n': 'n_opt_qn',
        'target': 'n_opt_qn_target'
    },
    'Tune set': {
        'n': 'n_opt_tune',
        'target': 'n_opt_tune_target'
    },
    'Scale set': {
        'n': 'n_opt_scale',
        'target': 'n_opt_scale_target'
    }
}
for prob_set in ['Basic set', 'Qn set', 'Tune set', 'Scale set']:

    # Subtitle
    print(prob_set)

    _counter_total = 0
    _counter_total_target = 0

    for prob in opt_problem_counters:

        _counter_total += opt_problem_counters[prob][set_info[prob_set]['n']]
        _counter_total_target += opt_problem_counters[prob][set_info[prob_set]['target']]

        # Compute percentage
        try:
            percent = opt_problem_counters[prob][set_info[prob_set]['n']] / \
                opt_problem_counters[prob][set_info[prob_set]['target']]*100.0
        except:
            percent = 0.0

        # Print stat
        print('prob: {:>26s}   {:4d}/{:4d}   ({:.2f}%)'.format(
            prob, opt_problem_counters[prob][set_info[prob_set]['n']],
            opt_problem_counters[prob][set_info[prob_set]['target']], percent))

    print('total:{:>26s}   {:4d}/{:4d}   ({:.2f}%)'.format('', _counter_total, _counter_total_target,
        _counter_total/_counter_total_target*100.0))

print('-------------------------------------------------------')
print('total:{:>26s}   {:4d}/{:4d}   ({:.2f}%)'.format('', counter_total, counter_total_target,
    counter_total/counter_total_target*100.0))