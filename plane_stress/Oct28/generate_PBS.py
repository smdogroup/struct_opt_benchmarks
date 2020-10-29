import argparse
import os

p = argparse.ArgumentParser()
p.add_argument('jobname', type=str, metavar='[jobname]')
p.add_argument('n', type=int, metavar='[number of PBS scripts]')
args = p.parse_args()

for i in range(args.n):
    PBSname=args.jobname + '_part{:d}.pbs'.format(i+1)
    with open(PBSname, 'w') as f:
        f.write('#PBS -N {:s}_part{:d}\n'.format(args.jobname, i+1))
        f.write('#PBS -l nodes=1:ppn=24\n')
        f.write('#PBS -l walltime=48:00:00\n')
        f.write('#PBS -l pmem=4gb\n')
        f.write('#PBS -q kennedy-lab\n')
        f.write('#PBS -o {:s}_part{:d}.out\n'.format(args.jobname, i+1))
        f.write('#PBS -j oe\n')
        f.write('\n')
        f.write('cd $PBS_O_WORKDIR\n')
        f.write('echo Start running job {:s}_part{:d}\n'.format(args.jobname, i+1))
        f.write('bash {:s}_part{:d}.sh\n'.format(args.jobname, i+1))
