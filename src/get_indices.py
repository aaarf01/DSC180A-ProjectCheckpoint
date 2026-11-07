"""Parse NOAA monthly climate indices into (decimal year, value) arrays.
ONI: https://psl.noaa.gov/data/correlation/oni.data   PDO: https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/index/ersst.v5.pdo.dat
usage: get_indices.py <oni.data> <pdo.dat> <data dir>"""
import sys, numpy as np
def monthly(path, skip):
    r = []
    for line in open(path).read().splitlines()[skip:]:
        p = line.split()
        if len(p) == 13 and p[0].isdigit():
            for mth, v in enumerate(p[1:], 1):
                v = float(v)
                if -90 < v < 90: r.append((int(p[0])+mth/12, v))
    return np.array(r)
oni, pdo, d = sys.argv[1:4]
np.save(f'{d}/oni.npy', monthly(oni, 1)); np.save(f'{d}/pdo.npy', monthly(pdo, 2))
