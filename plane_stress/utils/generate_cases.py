#!/usr/bin/env python3

"""
This script generates a batch of cases and an executable that submits all cases to PACE.
This script must be called in the directory which has a folder called pkls/, each pickle
file in pkls/ must have mesh, boundary condition and load defined

before calling this script:

main-folder
       |---pkls/

after calling this script:

main-folder
       |---pkls/
       |---cases.sh
       |---cases_part1.sh
       |---cases_part2.sh
            ...
       |---cases_partM.pbs
       |---cases_par1.pbs
       |---cases_part2.pbs
            ...
       |---cases_partM.pbs
       |---submit
       |---results/
               |---[case 1 placeholder folder]
               |---[case 2 placeholder folder]
                    ...
               |---[case n placeholder folder]
"""


import os
import re
import argparse
import plane_stress

# set up parser
p = argparse.ArgumentParser()
p.add_argument('npart', type=int, help='number of partitions')
p.add_argument('--opt_problem', nargs='*', type=str, default=None, choices=[
    'comp_min_mass_constr', 'comp_min_massfreq_constr',
    'comp_min_massstress_constr', 'comp_min_massfreqstress_constr',
    'stress_min_mass_constr', 'mass_min_stress_constr'])
p.add_argument('--optimizer', nargs='*', type=str, default=None, choices=[
    'all', 'ParOpt', 'ParOptAdapt', 'ParOptFilter',
    'ParOptFilterSoc', 'SNOPT', 'IPOPT'])
p.add_argument('--walltime', type=int, default=5, help='walltime requested in hours')
args = p.parse_args()

# Input value check
if args.opt_problem is None:
    raise ValueError("\n--opt_problem must be specified. See generate_cases.py --help for more info.")
if args.optimizer is None:
    raise ValueError("\n--optimizer must be specified. See generate_cases.py --help for more info.")

# Baseline frequencies
# Note that these are frequencies of final designs for mass constrained compliance
# minimization problems using SNOPT, mass < 0.4, n = 64, qval = 8.0
base_freqs = plane_stress.utils.base_freqs

# get all pkl file names
pkls = os.listdir('pkls')

# Optimization script
optimize = '~/git/struct_opt_benchmarks/plane_stress/optproblems/optimize.py'

# Set some problem parameters to be constant
design_mass = '--design_mass 0.4'
design_stress = '--design_stress 0.9 --stress_as_fraction'
qval_ks_epsilon_iter = '--qval 8.0 --epsilon 0.1 --ks_parameter 100.0 --max_iter 1000'

# Some other constants
freq_ratio = 1.2
total_required_cases = 0
existing_cases = 0
new_cases = 0

# Get optimization problem and optimizer list from input
opt_problems = args.opt_problem
if 'all' in args.optimizer:
    optimizers = ['ParOpt', 'ParOptAdapt', 'ParOptFilter',
        'ParOptFilterSoc', 'SNOPT', 'IPOPT']
else:
    optimizers = args.optimizer

# Define corresponding optimizer options here
optimizer_setting = {
    'ParOpt':'ParOpt',
    'ParOptAdapt':'ParOpt --ParOpt_use_adaptive_gamma_update',
    'ParOptFilter':'ParOpt --ParOpt_use_filter',
    'ParOptFilterSoc':'ParOpt --ParOpt_use_filter --ParOpt_use_soc',
    'SNOPT':'SNOPT',
    'IPOPT':'IPOPT'}

# Create the result folder
try:
    os.mkdir('results')
except:
    pass

# Create placeholder folders for all cases
for opt_problem in opt_problems:
    for pkl in pkls:

        # folder name = pkl name + opt_problem
        casefolder = 'results/{:s}-{:s}'.format(opt_problem, os.path.splitext(pkl)[0])

        # Make directory if does not exist
        try:
            os.mkdir(casefolder)
        except:
            pass

# get name of all existing files in results/ folder
casefolders = os.listdir('results')
current_files = []
for casefolder in casefolders:
    current_files.extend(os.listdir('results/'+casefolder))

# Generate a shell script containing all cases
with open('cases.sh', 'w') as f:
    for opt_problem in opt_problems:
        for optimizer in optimizers:
            for pkl in pkls:

                total_required_cases += 1

                # Write the case to runcases script only if there
                # is no result pkl file for that case

                caseExists = False
                for current_file in current_files:
                    if ('.pkl' in current_file and
                        opt_problem in current_file and
                        optimizer in current_file and
                        '-'+os.path.splitext(pkl)[0]+'-q' in current_file):
                        caseExists = True
                        existing_cases += 1
                        break

                if not caseExists:
                    # Get design frequency
                    meshtype = pkl.split('-')[0]
                    domain = pkl.split('-')[1]
                    AR = pkl.split('-')[3]
                    if 'hole' in pkl:
                        hole = 'hole'
                    else:
                        hole = 'nohole'

                    # set output directory
                    outdir = 'results/{:s}-{:s}'.format(opt_problem, os.path.splitext(pkl)[0])

                    if meshtype == 'structured':
                        design_freq = '--design_freq {:.3f}'.format(base_freqs[meshtype][domain][AR]*freq_ratio)
                    else:
                        design_freq = '--design_freq {:.3f}'.format(base_freqs[meshtype][domain][AR][hole]*freq_ratio)

                    if opt_problem == 'comp_min_mass_constr':
                        line = '{:s} pkls/{:s} {:s} {:s} --outdir {:s} {:s} {:s}\n'.format(
                            optimize, pkl, opt_problem, optimizer_setting[optimizer], outdir, qval_ks_epsilon_iter, design_mass)

                    elif opt_problem == 'comp_min_massfreq_constr':
                        line = '{:s} pkls/{:s} {:s} {:s} --outdir {:s} {:s} {:s} {:s}\n'.format(
                            optimize, pkl, opt_problem, optimizer_setting[optimizer], outdir, qval_ks_epsilon_iter, design_mass, design_freq)

                    elif opt_problem == 'comp_min_massstress_constr':
                        line = '{:s} pkls/{:s} {:s} {:s} --outdir {:s} {:s} {:s} {:s}\n'.format(
                            optimize, pkl, opt_problem, optimizer_setting[optimizer], outdir, qval_ks_epsilon_iter, design_mass, design_stress)

                    elif opt_problem == 'comp_min_massfreqstress_constr':
                        line = '{:s} pkls/{:s} {:s} {:s} --outdir {:s} {:s} {:s} {:s} {:s}\n'.format(
                            optimize, pkl, opt_problem, optimizer_setting[optimizer], outdir, qval_ks_epsilon_iter, design_mass, design_freq, design_stress)

                    elif opt_problem == 'stress_min_mass_constr':
                        line = '{:s} pkls/{:s} {:s} {:s} --outdir {:s} {:s} {:s} {:s}\n'.format(
                            optimize, pkl, opt_problem, optimizer_setting[optimizer], outdir, qval_ks_epsilon_iter, design_mass, design_stress)

                    elif opt_problem == 'mass_min_stress_constr':
                        line = '{:s} pkls/{:s} {:s} {:s} --outdir {:s} {:s} {:s}\n'.format(
                            optimize, pkl, opt_problem, optimizer_setting[optimizer], outdir, qval_ks_epsilon_iter, design_stress)

                    f.write(line)
                    new_cases += 1

# Partition run case script if specified:
if args.npart > 1:

    # read in all lines in cases.sh
    with open('cases.sh', 'r') as f:
        cases_sh_lines = f.readlines()

    # Generate partitioned scripts
    for i in range(args.npart):
        part_name = 'cases_part{:d}.sh'.format(i+1)

        # Create partitioned scripts
        with open(part_name, 'w') as f:
            f.write('#This is part {:d} / {:d}\n'.format(i+1, args.npart))

    # Write them cases_part{x}.sh
    counter = 0
    for line in cases_sh_lines:
        part_name = 'cases_part{:d}.sh'.format(counter % args.npart + 1)

        with open(part_name, 'a') as f:
            f.write(line)

        counter += 1

    # Generate PBS script
    for i in range(args.npart):
        pbs_name = 'cases_part{:d}.pbs'.format(i+1)

        with open(pbs_name, 'w') as f:
            f.write('#PBS -N part{:d}\n'.format(i+1))
            f.write('#PBS -A GT-gkennedy9-CODA20\n')
            f.write('#PBS -l nodes=1:ppn=24\n')
            f.write('#PBS -l walltime={:d}:00:00\n'.format(args.walltime))
            f.write('#PBS -l pmem=4gb\n')
            f.write('#PBS -o part{:d}.out\n'.format(i+1))
            f.write('#PBS -j oe\n')
            f.write('#PBS -m abe\n')
            f.write('#PBS -M yfu97@gatech.edu\n')
            f.write('\n')
            f.write('cd $PBS_O_WORKDIR\n')
            f.write('echo Job starts at:\n')
            f.write('date\n')
            f.write('parallel -j 24 :::: cases_part{:d}.sh\n'.format(i+1))
            f.write('echo Job ends at:\n')
            f.write('date\n')
            f.write('date >> done_part{:d}\n'.format(i+1))

# Generate a script that submits all PBS scripts
with open('submit', 'w') as f:
    for i in range(args.npart):
        f.write('qsub cases_part{:d}.pbs\n'.format(i+1))

# Make it executable
os.system('chmod +x submit')

# Print summary:
print("-------------------------------")
print("Done generating cases, summary:")
print("Total number of cases   : {:d}".format(total_required_cases))
print("Number of existing cases: {:d}".format(existing_cases))
print("Number of new cases     : {:d}".format(new_cases))
print("-------------------------------")
