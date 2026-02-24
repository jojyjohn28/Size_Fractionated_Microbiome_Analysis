#!/usr/bin/env python3
import os, sys
import pandas as pd

idx_dir = sys.argv[1]          # idxstats folder
samples_fp = sys.argv[2]       # samples.txt (one sample per line, without .idxstats.tsv)
out_fp = sys.argv[3]           # output tsv

samples = [x.strip() for x in open(samples_fp) if x.strip()]
df = None

for i, s in enumerate(samples, 1):
    fp = os.path.join(idx_dir, f"{s}.idxstats.tsv")
    if not os.path.exists(fp):
        print(f"[WARN] missing: {fp}", file=sys.stderr)
        continue

    # idxstats columns: ref, length, mapped, unmapped
    t = pd.read_csv(fp, sep="\t", header=None, names=["id","len","mapped","unmapped"])
    t = t[t["id"] != "*"][["id","mapped"]].rename(columns={"mapped": s})

    df = t if df is None else df.merge(t, on="id", how="outer")

    if i % 5 == 0:
        print(f"[INFO] processed {i}/{len(samples)}", file=sys.stderr)

df = df.fillna(0)

# keep integers
for c in df.columns[1:]:
    df[c] = df[c].astype("int64")

df.to_csv(out_fp, sep="\t", index=False)
print(f"[DONE] wrote {out_fp} with {df.shape[0]} genes and {df.shape[1]-1} samples", file=sys.stderr)

