"""Run once woa23_all_n00_01.nc is available (NOAA server was down). NOT yet run on real data.
Gives Fig 1A (paper r=-0.52), Fig 2B field, and Fig 4D with a real nutricline. usage: woa_steps.py ../data ../out"""
import sys, numpy as np, netCDF4 as nc, lib
d, o = sys.argv[1:3]
w = nc.Dataset(f'{d}/woa23_all_n00_01.nc'); N = np.ma.filled(w['n_an'][0], np.nan)*1.025e-3*1e3  # umol/kg to uM, density 1.025 assumed
Zd = w['depth'][:]; lat = w['lat'][:]; lon = w['lon'][:]
def zno3(p, thr=3.):
    ok = np.isfinite(p)
    if ok.sum() < 2: return np.nan
    x, y = Zd[ok], p[ok]
    if y[0] >= thr: return 0.
    i = np.argmax(y >= thr)
    return np.nan if y[i] < thr else x[i-1]+(thr-y[i-1])*(x[i]-x[i-1])/(y[i]-y[i-1])
Z = np.apply_along_axis(zno3, 0, N)                        # lat S->N, lon -180->180
Z = np.concatenate([Z[:, 200:], Z[:, :200]], 1)[::-1][49:131].ravel('F')  # match theta_LL layout
s = np.load(f'{o}/sat_mask.npz'); m = s['m']; zz = Z[m]; k = np.isfinite(zz)
print('Fig 1A r', np.corrcoef(zz[k], s['mean'][k])[0, 1], 'n', k.sum())
np.save(f'{o}/woa_zno3.npy', Z); T = s['T']; la = s['LAT'][m]
for n, q in (('All', k), ('All lowNut', k & (zz > 50)), ('All highNut', k & (zz < 50)), ('N lowNut', k & (zz > 50) & (la > 0)), ('N highNut', k & (zz < 50) & (la > 0)),
             ('S lowNut', k & (zz > 50) & (la < 0)), ('S highNut', k & (zz < 50) & (la < 0))):
    print(n, q.sum(), lib.boot_ci_mean(T[q]))
