import sys, numpy as np
from deeptools import heatmapper
print("sample\tprof_max\tprof_min\theat_p2\theat_p98")
for f in sys.argv[1:]:
    hm = heatmapper.heatmapper(); hm.read_matrix_file(f)
    M = hm.matrix.matrix; b = hm.matrix.sample_boundaries
    for i, l in enumerate(hm.matrix.sample_labels):
        x = M[:, b[i]:b[i+1]]; p = np.nanmean(x, axis=0)
        print(f"{l}\t{p.max():.3f}\t{p.min():.3f}\t{np.nanpercentile(x,2):.3f}\t{np.nanpercentile(x,98):.3f}")
