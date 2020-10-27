from plane_stress_analysis import PlaneStressAnalysis
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import argparse
import pickle

# Set up parser
p = argparse.ArgumentParser()
p.add_argument('filename', metavar='cantilever.pkl', type=str)
args = p.parse_args()

# Load in pickle file
with open(args.filename, 'rb') as pklfile:
    prob_pkl = pickle.load(pklfile)

# Get data
prob_name = prob_pkl['prob_name']
nelems = prob_pkl['nelems']
nnodes = prob_pkl['nnodes']
nvars = prob_pkl['nvars']
C = prob_pkl['C']
conn = prob_pkl['conn']
X = prob_pkl['X']
vars = prob_pkl['vars']
force = prob_pkl['force']
r0 = prob_pkl['r0']
x = prob_pkl['x']

# Check if we have a solution in pickle or not
if x is None:
    np.random.seed(0)
    x = np.random.rand(nnodes)

# Instantiate analysis because we need to compute filtered design
qval = 5.0
density = 2700.0
freq= 0.6
lambda0 = (2.0*np.pi*freq)**2
ks = 100.0
num_eigs = 8
sigma = -100.0
analysis = PlaneStressAnalysis(conn, vars, X, force, r0, qval, C,
    density=density, lambda0=lambda0, ks_parameter=ks, num_eigs=num_eigs,
    eigshsigma=sigma)

analysis.plot_solution(x)