import argparse
import os

p = argparse.ArgumentParser()
p.add_argument('jobsfile', type=str, help='a single file contains all jobs')
p.add_argument('n', type=int, help='number of files to split')
args = p.parse_args()

with open(args.jobsfile, 'r') as f:
    lines = f.readlines()

njobs = len(lines)
nfiles = args.n

njob = [ njobs // nfiles ] * nfiles

for i in range(njobs % nfiles):
    njob[i] += 1

index = 1
start = 0
for i in range(nfiles):
    filename = os.path.splitext(args.jobsfile)[0] + \
        '_part{:d}'.format(index) + os.path.splitext(args.jobsfile)[1]
    with open(filename, 'w') as f:
        for line in range(njob[i]):
            text = lines[start+line]
            oldoutdir = text.split()[text.split().index('--outdir')+1]
            newoutdir = oldoutdir + '_part{:d}'.format(index)
            text = text.replace('--outdir {:s}'.format(oldoutdir),'--outdir {:s}'.format(newoutdir))
            f.write(text)
        start += njob[i]
    index += 1
