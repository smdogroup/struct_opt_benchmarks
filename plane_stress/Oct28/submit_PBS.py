import argparse

p = argparse.ArgumentParser()
p.add_argument('jobname', type=str, metavar='[jobname]')
p.add_argument('n', type=int, metavar='[number of PBS scripts]')
args = p.parse_args()

with open('submit_PBS.sh', 'w') as f:
    for i in range(args.n):
        f.write('qsub {:}_part{:d}.pbs\n'.format(args.jobname, i+1))
