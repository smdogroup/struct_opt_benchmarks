fname = "generate_pkls.sh"

lxs = [1.0, 2.0, 3.0]
nys = [40, 80, 120]

with open(fname, 'w') as f:
    for i in range(len(lxs)):
        for j in range(len(nys)):
            lx = lxs[i]
            ly = 1.0
            ny = nys[j]
            nx = int(ny * lx / ly)
            line = "python3 ../preproc-cantilever.py --outdir pkls --nx {:d} --ny {:d} --lx {:.1f} --ly {:.1f}\n".format(
                nx, ny, lx, ly)
            f.write(line)