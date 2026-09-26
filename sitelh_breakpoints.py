#!/usr/bin/env python3
"""
Locate a recombinant block from IQ-TREE per-site log-likelihoods under two fixed trees.
 
Usage:
  python3 sitelh_breakpoints.py vert.sitelh rec.sitelh [--perm 9999] [--plot out.png] [--offset 0]
  python3 sitelh_breakpoints.py vert.sitelh rec.sitelh --perm 0 --tsv     # one line: stat start end
 
Inputs are the .sitelh files written by
  iqtree3 -s full.fasta -m MODEL -te vert.treefile --sitelh --prefix vert
  iqtree3 -s full.fasta -m MODEL -te rec.treefile  --sitelh --prefix rec
 
d[i] = lnL_i(recombinant tree) - lnL_i(vertical tree). The recombinant block is the
contiguous run of columns [a, b] maximising  sum(d inside) - sum(d outside), found
exactly in O(n) from the cumulative sum. The statistic is that maximum.
 
Significance, two options:
  --perm N   site-permutation test (quick; asks whether sites favouring the recombinant
             tree are more clustered than chance)
  parametric bootstrap (stricter; see run_parametric_bootstrap.sh): run this script
             with --perm 0 --tsv on alignments simulated under the vertical tree.
"""
import argparse
import sys
import numpy as np
 
 
def read_sitelh(path):
    lines = [l for l in open(path).read().split("\n") if l.strip()]
    if len(lines) < 2 or not lines[1].startswith("Site_Lh"):
        sys.exit(f"{path}: not a single-tree .sitelh file")
    return np.array([float(x) for x in lines[1].split()[1:]])
 
 
def best_block(d):
    """Return (statistic, start, end); start/end 1-based inclusive."""
    c = np.concatenate([[0.0], np.cumsum(d)])
    best, bi, bj, mn, mi = -np.inf, 0, 0, 0.0, 0
    for j in range(1, len(c)):
        if c[j] - mn > best:
            best, bi, bj = c[j] - mn, mi, j
        if c[j] < mn:
            mn, mi = c[j], j
    return 2 * best - c[-1], bi + 1, bj
 
 
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vertical")
    ap.add_argument("recombinant")
    ap.add_argument("--perm", type=int, default=9999)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--offset", type=int, default=0,
                    help="add to reported coordinates if the alignment starts after column 1")
    ap.add_argument("--plot")
    ap.add_argument("--tsv", action="store_true", help="print only: statistic start end")
    a = ap.parse_args()
 
    v, r = read_sitelh(a.vertical), read_sitelh(a.recombinant)
    if len(v) != len(r):
        sys.exit("the two .sitelh files have different numbers of sites")
    d = r - v
    n = len(d)
    stat, s, e = best_block(d)
    if a.tsv:
        print(f"{stat:.4f}\t{s + a.offset}\t{e + a.offset}")
        return
 
    print(f"sites: {n}")
    print(f"recombinant block: {s + a.offset}-{e + a.offset}")
    print(f"sum dlnL inside = {d[s-1:e].sum():.2f}; outside = {d[:s-1].sum() + d[e:].sum():.2f}; "
          f"whole alignment = {d.sum():.2f}")
    print(f"block statistic = {stat:.2f}")
    if a.perm > 0:
        rng = np.random.default_rng(a.seed)
        null = np.array([best_block(rng.permutation(d))[0] for _ in range(a.perm)])
        p = (1 + np.sum(null >= stat)) / (a.perm + 1)
        print(f"site-permutation p = {p:.4g} ({a.perm} permutations)")
    else:
        p = None
 
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        x = np.arange(1, n + 1) + a.offset
        fig, ax = plt.subplots(figsize=(8, 3.2))
        ax.plot(x, np.cumsum(d), color="#2a6f97", lw=1.5)
        ax.axvspan(s + a.offset, e + a.offset, color="#e9c46a", alpha=0.35, lw=0)
        ax.axhline(0, color="#888", lw=0.6)
        ax.set_xlabel("alignment column")
        ax.set_ylabel("cumulative ΔlnL\n(recombinant − vertical)")
        title = f"Inferred block {s + a.offset}–{e + a.offset}"
        if p is not None:
            title += f", permutation p = {p:.3g}"
        ax.set_title(title, fontsize=10)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        fig.tight_layout()
        fig.savefig(a.plot, dpi=200)
        print(f"plot written to {a.plot}")
 
 
if __name__ == "__main__":
    main()