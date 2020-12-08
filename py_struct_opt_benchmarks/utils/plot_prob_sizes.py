#!/usr/bin/env python3

import os
import argparse
import pickle
from pprint import pprint
import matplotlib.pyplot as plt
import numpy as np

p = argparse.ArgumentParser()
p.add_argument('--pkl_folder', default='pkls')
args = p.parse_args()

# Get pkls
pkls = os.listdir(args.pkl_folder)

# Get domains
domains = []

for pkl in pkls:
    if pkl.split('-')[0] == '3D':
        domain = pkl.split('-')[2]
    else:
        domain = pkl.split('-')[1]

    if domain not in domains:
        domains.append(domain)

domains.sort()

if len(domains) == 4:
    domains = ['cantilever', 'michell', 'MBB', 'lbracket']


# Re-arrange pkls list
pkls_reordered = []
for domain in domains:
    _pkls = []
    for pkl in pkls:
        if domain in pkl:
            _pkls.append(pkl)
    _pkls.sort()
    pkls_reordered += _pkls

# Loop over to get dof for each pickle
dof = []
dividers = []
domains = []
index = 0
for pkl in pkls_reordered:

    # Check if problem is 3D
    if pkl.split('-')[0] == '3D':
        dimension = 3
        domain = pkl.split('-')[2]
    else:
        dimension = 2
        domain = pkl.split('-')[1]

    if domain not in domains:
        domains.append(domain)
        dividers.append(index)

    with open('{:s}/{:s}'.format(args.pkl_folder, pkl), 'rb') as f:
        dof.append(pickle.load(f)['nnodes']*dimension)

    index += 1
dividers.append(index)

# Set up plotting environment
mpl_style_path = os.path.dirname(os.path.realpath(__file__)) + '/paper.mplstyle'
plt.style.use(mpl_style_path)
fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8))

# Rename legends
legends = {
    'cantilever': 'cantilever',
    'michell': 'Michell',
    'MBB': 'MBB',
    'lbracket': 'L-bracket'
}

# Plot
indices = list(range(len(dof)))

dof = [entry/1000 for entry in dof]
for i in range(len(domains)):
    d1 = dividers[i]
    d2 = dividers[i+1]
    ax.bar(indices[d1:d2], dof[d1:d2], label=legends[domains[i]])

ax.set_xlabel('Problem index')
ax.set_ylabel(r'Degrees of freedom, $\times 10^3$')

plt.legend()

if dimension == 3:
    plt.savefig('3D-prob_sizes.pdf')
    print('3D max dof: {:d}'.format(int(max(dof)*1000)))
    print('3D min dof: {:d}'.format(int(min(dof)*1000)))
else:
    plt.savefig('2D-prob_sizes.pdf')
    print('2D max dof: {:d}'.format(int(max(dof)*1000)))
    print('2D min dof: {:d}'.format(int(min(dof)*1000)))
plt.close()