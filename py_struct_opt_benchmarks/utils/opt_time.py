#!/usr/bin/env python3

import argparse
import pickle

def opt_time(pkl):

    with open(pkl, 'rb') as pklfile:
        prob_pkl = pickle.load(pklfile)

    print(prob_pkl['opt_time'])
    return prob_pkl['opt_time']

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('pkl', metavar=['pkl file'], type=str)
    args = p.parse_args()
    opt_time(args.pkl)


