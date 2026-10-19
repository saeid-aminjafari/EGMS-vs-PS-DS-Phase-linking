#!/usr/bin/env python3
"""
flag_tep_adi_points.py

Identify SARvey points that entered the solution ONLY through the MiaplPy PS
mask (TEP + amplitude dispersion), i.e. points inside maskPS.h5 whose temporal
phase coherence (TPC) is below the threshold. These are the points that would
not exist if SARvey had been run with use_ps: false.

Also checks whether the PS mask affected the FIRST-order network (p1). If it did
not, removing the flagged second-order points afterwards is equivalent to
reprocessing, because SARvey densifies each second-order point individually
against the first-order network.

Usage (run once per corridor):
    python3 flag_tep_adi_points.py \
        --p2   outputs/p2_coh70_ts.h5 \
        --p1   outputs/p1_ts.h5 \
        --tcoh outputs/temporal_coherence.h5 \
        --mask inverted/maskPS.h5 \
        --thr  0.7 --thr-p1 0.9 \
        --out  flagged_points.csv

Dataset names inside the .h5 files are detected automatically; use --list to
print what each file contains if detection fails.
"""
import argparse
import csv
import sys

import h5py
import numpy as np


def list_file(path):
    print("\n%s" % path)
    with h5py.File(path, "r") as f:
        def show(name, obj):
            if isinstance(obj, h5py.Dataset):
                print("   %-30s %-20s %s" % (name, str(obj.shape), obj.dtype))
        f.visititems(show)


def read_2d(path, preferred=None):
    """Return the 2-D raster in a file (mask or coherence)."""
    with h5py.File(path, "r") as f:
        if preferred and preferred in f:
            return np.asarray(f[preferred][()])
        cands = []
        f.visititems(lambda n, o: cands.append(n)
                     if isinstance(o, h5py.Dataset) and o.ndim == 2 else None)
        if not cands:
            sys.exit("No 2-D dataset found in %s (use --list)" % path)
        if len(cands) > 1:
            print("   note: %s has several 2-D datasets %s; using '%s'"
                  % (path, cands, cands[0]))
        return np.asarray(f[cands[0]][()])


def read_points(path):
    """Return (coord_xy [n,2] int, point_id [n] or None, lalo [n,2] or None)."""
    with h5py.File(path, "r") as f:
        if "coord_xy" not in f:
            sys.exit("'coord_xy' not found in %s (use --list)" % path)
        xy = np.asarray(f["coord_xy"][()]).astype(int)
        pid = np.asarray(f["point_id"][()]) if "point_id" in f else None
        lalo = np.asarray(f["coord_lalo"][()]) if "coord_lalo" in f else None
    return xy, pid, lalo


def lookup(raster, xy, order):
    a, b = (xy[:, 0], xy[:, 1]) if order == "rc" else (xy[:, 1], xy[:, 0])
    ok = (a >= 0) & (a < raster.shape[0]) & (b >= 0) & (b < raster.shape[1])
    out = np.full(len(xy), np.nan)
    out[ok] = raster[a[ok], b[ok]]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p2", required=True)
    ap.add_argument("--p1", required=True)
    ap.add_argument("--tcoh", required=True)
    ap.add_argument("--mask", required=True)
    ap.add_argument("--thr", type=float, default=0.7)
    ap.add_argument("--thr-p1", type=float, default=0.9)
    ap.add_argument("--out", default="flagged_points.csv")
    ap.add_argument("--list", action="store_true", help="print file contents and exit")
    a = ap.parse_args()

    if a.list:
        for p in (a.p2, a.p1, a.tcoh, a.mask):
            list_file(p)
        return

    tcoh = read_2d(a.tcoh).astype(float)
    mask = read_2d(a.mask) > 0
    if tcoh.shape != mask.shape:
        sys.exit("Shape mismatch: temporal coherence %s vs maskPS %s"
                 % (tcoh.shape, mask.shape))

    p1_xy, _, _ = read_points(a.p1)
    p2_xy, p2_id, p2_lalo = read_points(a.p2)

    # Determine index order (row/col vs col/row) using p1: all p1 must pass thr_p1
    best = None
    for order in ("rc", "cr"):
        frac = np.nanmean(lookup(tcoh, p1_xy, order) >= a.thr_p1)
        if best is None or frac > best[1]:
            best = (order, frac)
    order, frac = best
    print("coord_xy interpreted as %s; %.1f%% of p1 points have TPC >= %.2f"
          % ("(row, col)" if order == "rc" else "(col, row)", 100 * frac, a.thr_p1))
    if frac < 0.95:
        print("   WARNING: fewer p1 points pass than expected. Check that the"
              " temporal coherence file is the one SARvey used.")

    # First-order network: did maskPS add points below thr_p1?
    p1_t = lookup(tcoh, p1_xy, order)
    p1_m = lookup(mask.astype(float), p1_xy, order) > 0
    p1_bad = p1_m & (p1_t < a.thr_p1)
    print("\nFIRST-ORDER (p1): %d points, %d inside maskPS, %d inside maskPS with"
          " TPC < %.2f" % (len(p1_xy), p1_m.sum(), p1_bad.sum(), a.thr_p1))
    if p1_bad.sum() == 0:
        print("   -> PS mask did not alter the first-order network. Post-hoc removal"
              " of flagged p2 points is equivalent to reprocessing.")
    else:
        print("   -> PS mask DID add first-order points. Post-hoc removal is NOT"
              " strictly equivalent; reprocessing from step 0 is advisable.")

    # Second-order points entering only through the PS mask
    p2_t = lookup(tcoh, p2_xy, order)
    p2_m = lookup(mask.astype(float), p2_xy, order) > 0
    flagged = p2_m & (p2_t < a.thr)
    n = len(p2_xy)
    print("\nSECOND-ORDER (p2): %d points" % n)
    print("   inside maskPS                        : %6d (%.1f%%)" % (p2_m.sum(), 100 * p2_m.mean()))
    print("   inside maskPS AND TPC < %.2f (flag)   : %6d (%.1f%%)"
          % (a.thr, flagged.sum(), 100 * flagged.mean()))
    print("   kept after removal                   : %6d" % (n - flagged.sum()))

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        hdr = ["point_id", "tpc", "in_maskPS", "flag_remove"]
        if p2_lalo is not None:
            hdr += ["lat", "lon"]
        w.writerow(hdr)
        for i in range(n):
            row = [int(p2_id[i]) if p2_id is not None else i,
                   "%.4f" % p2_t[i], int(p2_m[i]), int(flagged[i])]
            if p2_lalo is not None:
                row += ["%.7f" % p2_lalo[i, 0], "%.7f" % p2_lalo[i, 1]]
            w.writerow(row)
    print("\nWritten: %s  (join on point_id in QGIS, drop flag_remove = 1)" % a.out)


if __name__ == "__main__":
    main()