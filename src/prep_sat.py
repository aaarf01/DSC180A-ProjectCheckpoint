"""Faithful port of preprocess_LowLat_thetaPrime_v3_bathy.m + prepare_time_files.m."""
import sys, numpy as np, pandas as pd, lib
d = sys.argv[1]
B = np.load(f'{d}/theta.npy').reshape(180, 360, 897, order='F')   # lat N->S, lon -180->180
yrdy = np.load(f'{d}/theta_time.npy'); yr, yd = yrdy
dates = pd.to_datetime([f'{y}{x:03d}' for y, x in zip(yr, yd)], format='%Y%j')
timenum = (dates-pd.Timestamp('2000-01-01')).days.values.astype(float)
np.savez(f'{d}/time.npz', yr=yr, yd=yd, month=dates.month.values, timenum=timenum, time=yr+yd/365)
lat = 89.5-np.arange(180.); lon = np.arange(-179.5, 180)
def shift(A): return np.concatenate([A[:, 200:], A[:, :200]], 1)   # lon 20.5..379.5
B = shift(B)[49:131]
print('band', B.shape)
L = B.reshape(29520, 897, order='F')
print('cells with any data', np.isfinite(L).any(1).sum())
print('cells passing NaN filter (<=224 NaN)', (np.isnan(L).sum(1) <= 224).sum())
S = lib.preprocess(L, np.ones(29520, bool))
print('kept', np.isfinite(S).all(1).sum(), 'min/max', np.nanmin(S), np.nanmax(S))
np.save(f'{d}/theta_LL.npy', S.astype(np.float32))
