import sys, gzip, json, numpy as np
LOW = 10
seen = {}; nsamp = 0
for f in sys.argv[1:]:
    fh = gzip.open(f, "rt"); h = json.loads(fh.readline()[1:])
    labs, bnd = h["sample_labels"], h["sample_boundaries"]
    bs, up = h["bin size"][0], h["upstream"][0]; nsamp += len(labs)
    for line in fh:
        p = line.rstrip("\n").split("\t")
        v = np.nan_to_num(np.array(p[6:], dtype=float), nan=0.0)
        tss, strand = int(p[1]), p[5]
        for i, l in enumerate(labs):
            x = np.abs(v[bnd[i]:bnd[i+1]]); k = np.where(x > LOW)[0]
            if len(k) == 0: continue
            o1, o2 = -up + k.min() * bs, -up + (k.max() + 1) * bs
            s, e = (tss + o1, tss + o2) if strand == "+" else (tss - o2, tss - o1)
            r = seen.setdefault(p[3], [p[0], max(s, 0), e, strand, []])
            r[1], r[2] = min(r[1], max(s, 0)), max(r[2], e); r[4].append((l, float(x.max())))
print(f"{nsamp} samples scanned; regions with any 10-bp bin |value| above:")
for t in (10, 20, 50, 100):
    print(f"  > {t}: {sum(1 for r in seen.values() if max(m for _, m in r[4]) > t)} regions")
print("\ngene\tregion\tstrand\tn_samples>10\tmax\ttop samples")
rows = sorted(seen.items(), key=lambda kv: (-len(kv[1][4]), -max(m for _, m in kv[1][4])))
with open("flagged_regions_gt10.bed", "w") as b:
    for g, (c, s, e, st, hits) in rows:
        top = sorted(hits, key=lambda t: -t[1])[:3]
        print(f"{g}\t{c}:{s}-{e}\t{st}\t{len(hits)}\t{max(m for _, m in hits):.1f}\t" + " ".join(f"{l}={m:.0f}" for l, m in top))
        b.write(f"{c}\t{s}\t{e}\t{g}\t{len(hits)}\t{st}\n")
