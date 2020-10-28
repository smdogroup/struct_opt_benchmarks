from plane_stress_analysis import PlaneStressAnalysis
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import argparse
import pickle

# Set up parser
p = argparse.ArgumentParser()
p.add_argument('filename', metavar='[picklefile]', type=str)
args = p.parse_args()

# Load in pickle file
with open(args.filename, 'rb') as pklfile:
    prob_pkl = pickle.load(pklfile)

# Get data
prob_name = prob_pkl['prob_name']
nelems = prob_pkl['nelems']
nnodes = prob_pkl['nnodes']
ndof = prob_pkl['ndof']
C = prob_pkl['C']
conn = prob_pkl['conn']
X = prob_pkl['X']
dof = prob_pkl['dof']
force = prob_pkl['force']
r0 = prob_pkl['r0']
density = prob_pkl['density']
qval = prob_pkl['qval']
x = prob_pkl['x']
opt_settings = prob_pkl['opt_settings']

# Check if we have a solution in pickle or not
if x is None:
    np.random.seed(0)
    x = np.random.rand(nnodes)

# Instantiate analysis because we need to compute filtered design
if opt_settings is None:
    print("[Warning] pkl file doesn't contain optimization settings, use default values")
    analysis = PlaneStressAnalysis(conn, dof, X, force, r0, qval, C)
else:
    freq= opt_settings['freq']
    ks = opt_settings['eig_ks']
    num_eigs = opt_settings['num_eigs']
    sigma = opt_settings['eigshsigma']
    lambda0 = (2.0*np.pi*freq)**2
    analysis = PlaneStressAnalysis(conn, dof, X, force, r0, qval, C,
        density=density, lambda0=lambda0, ks_parameter=ks, num_eigs=num_eigs,
        eigshsigma=sigma)

analysis.plot_solution(x)