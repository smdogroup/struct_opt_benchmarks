import os

fname = "freq_cases.sh"

pkls = os.listdir('pkls')
methods = [
    '--optimizer ParOpt',
    '--optimizer ParOpt --ParOpt_use_adaptive_gamma_update',
    '--optimizer ParOpt --ParOpt_use_filter',
    '--optimizer ParOpt --ParOpt_use_filter --ParOpt_use_soc',
    '--optimizer IPOPT',
    '--optimizer SNOPT',
]
base_freqs = [0.382, 0.176, 0.0959]
masses = [0.3, 0.4, 0.5]
ratios = [1.2, 1.4, 1.6]

with open(fname, 'w') as f:
    for j in range(len(methods)):
        for k in range(len(masses)):
            for m in range(len(ratios)):
                for i in range(len(pkls)):
                    pkl = pkls[i]
                    if 'lx1.0' in pkl:
                        freq = base_freqs[0]*ratios[m]
                    elif 'lx2.0' in pkl:
                        freq = base_freqs[1]*ratios[m]
                    elif 'lx3.0' in pkl:
                        freq = base_freqs[2]*ratios[m]
                    line = 'python3 ../comp_min_massfreq_constr.py pkls/{:s} {:s} --mass {:.3f} --freq {:.3f} --max_iter 600 --outdir comp_min_massfreq_constr\n'.format(
                        pkl, methods[j], masses[k], freq)
                    f.write(line)
