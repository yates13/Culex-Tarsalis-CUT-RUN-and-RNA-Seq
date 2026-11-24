import sys, gzip, numpy as np
from deeptools import heatmapper
f, site = sys.argv[1], sys.argv[2]
hm = heatmapper.heatmapper(); hm.read_matrix_file(f)
M = hm.matrix.matrix; b = hm.matrix.sample_boundaries; labs = hm.matrix.sample_labels
rows = [l.split("\t")[:6] for l in gzip.open(f, "rt") if not l.startswith("@")]
print(f"{len(rows)} regions, {M.shape[1]} columns, samples: {len(labs)}")
EDGE = 10
for i, l in enumerate(labs):
    if not l.startswith(site): continue
    x = M[:, b[i]:b[i+1]]; e = np.nanmean(x[:, -EDGE:], axis=1)
    o = np.argsort(-np.nan_to_num(e, nan=-1e9))
    print(f"\n## {l}: edge mean all={np.nanmean(e):.3f}  minus top1={np.nanmean(np.delete(e, o[:1])):.3f}  minus top5={np.nanmean(np.delete(e, o[:5])):.3f}  TSS-bin mean={np.nanmean(x[:, x.shape[1]//2]):.3f}")
    for r in o[:5]:
        print(f"   {rows[r][0]}:{rows[r][1]}-{rows[r][2]} {rows[r][3]} strand={rows[r][5]} edge={e[r]:.1f}")
    top = o[0]
    print("   top region's edge value in every sample: " + "  ".join(f"{labs[j]}={np.nanmean(M[top, b[j]:b[j+1]][-EDGE:]):.1f}" for j in range(len(labs))))
