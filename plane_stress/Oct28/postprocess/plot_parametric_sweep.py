import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
import numpy as np
import argparse
import glob
import ntpath
import os

p = argparse.ArgumentParser()
p.add_argument('datafolder', type=str)
p.add_argument('sweep', type=str, choices=[
    'freq', 'mass', 'AR', 'meshsize'])
p.add_argument('--freqratio', type=float, default=1.4, choices=[1.2, 1.4, 1.6])
p.add_argument('--mass', type=float, default=0.4, choices=[0.3, 0.4, 0.5])
p.add_argument('--AR', type=int, default=2, choices=[1, 2, 3])
p.add_argument('--ny', type=int, default=60, choices=[30, 60, 90])
p.add_argument('--save', action='store_true')
args = p.parse_args()

base_freqs = [0.382, 0.176, 0.0959]
ratios = [1.2, 1.4, 1.6]
masses = [0.3, 0.4, 0.5]
ARs = [1, 2, 3]
nys = [30, 60, 90]
methods = ['ParOpt', 'ParOptadapt', 'ParOptfilter', 'ParOptfiltersoc', 'IPOPT', 'SNOPT']

pngs = [ ntpath.basename(name) for name in glob.glob("{:s}/*.png".format(args.datafolder))]

# Create a figure
plt.figure(figsize=(6,8), dpi=300)

# Create a grid
grid = gridspec.GridSpec(nrows=7, ncols=7)
grid.update(wspace=0, hspace=0)

mass = [args.mass]*3
freq = [args.freqratio*base_freqs[ARs.index(args.AR)]]*3
nx = [args.AR*args.ny]*3
ny = [args.ny]*3
ly = 1.0
lx = [ly*args.AR]*3
AR = [args.AR]*3

if args.sweep == 'freq':
    freq = [r*base_freqs[ARs.index(args.AR)] for r in ratios]
elif args.sweep == 'mass':
    mass = masses
elif args.sweep == 'AR':
    AR = ARs
    nx = [args.ny*AR for AR in ARs]
    lx = [ly*AR for AR in ARs]
elif args.sweep == 'meshsize':
    ny = nys

for row in range(3):
    for col in range(6):
        fname = 'comp_min_massfreq_constr-{:s}-mass-{:.3f}-freq-{:.3f}-cantilever-nx{:d}-ny{:d}-lx{:.1f}-ly{:.1f}'.format(
                methods[col], mass[row], freq[row], nx[row], ny[row], lx[row], ly)

        if fname+'.png' in pngs:
            ax = plt.subplot(grid[2*row+1, col+1])
            ax.axis('off')
            img = mpimg.imread('{:s}/{:s}.png'.format(args.datafolder, fname))
            ax.imshow(img)

        if fname+'-history.png' in pngs:
            ax = plt.subplot(grid[2*row+2, col+1])
            ax.axis('off')
            img = mpimg.imread('{:s}/{:s}-history.png'.format(args.datafolder, fname))
            ax.imshow(img)

# Print table head
for col in range(6):
    ax = plt.subplot(grid[0, col+1])
    ax.axis('off')
    ax.text(0,0,methods[col], fontsize=5)

for row in range(3):
    ax = plt.subplot(grid[2*row+1, 0])
    ax.axis('off')
    ax.text(0,0.3, 'mass  = {:.3f}'.format(mass[row]),fontsize=5)
    ax.text(0,0.1, 'freq  = {:.3f}'.format(freq[row]),fontsize=5)
    ax.text(0,-0.1,'AR    = {:d}'.format(AR[row]),fontsize=5)
    ax.text(0,-0.3,'ny    = {:d}'.format(ny[row]),fontsize=5)

if (args.save):
    plt.savefig('sweep-{:s}-mass-{:.3f}-freq-{:.3f}-AR-{:d}-ny{:d}.png'.format(
        args.sweep, args.mass, args.freq, args.AR, args.ny))
else:
    plt.show()
