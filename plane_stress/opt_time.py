import argparse
import pickle

p = argparse.ArgumentParser()
p.add_argument('pkl', metavar=['pkl file'], type=str)
args = p.parse_args()

with open(args.pkl, 'rb') as pklfile:
    prob_pkl = pickle.load(pklfile)

print(prob_pkl['opt_time'])

