"""Figs 2-4, S7-S16 numbers. Ports ANOVA_*.m, EOF_variance.m, longterm_SST_theta_overlay.m, delta_SST_theta_plot.m."""
import sys, json, numpy as np, pandas as pd, lib
d, o = sys.argv[1:3]
t = np.load(f'{d}/time.npz'); yr, month, timenum, time = t['yr'], t['month'], t['timenum'], t['time']
yr20 = np.where(yr == 2022, 2021, yr)       # one composite dated 2022-001 joined to 2021 so years=20 (LOG A10)
TH = np.load(f'{d}/theta_LL.npy').astype(float)
SS = np.load(f'{d}/SST.npy').reshape(180, 360, 897, order='F')
SS = np.concatenate([SS[:, 200:], SS[:, :200]], 1)[49:131].reshape(29520, 897, order='F')
la = (40.5-np.arange(82.))[:, None]*np.ones((1, 360)); lo = np.ones((82, 1))*(20.5+np.arange(360.))[None]
LAT, LON = la.ravel('F'), lo.ravel('F'); W = np.cos(np.radians(LAT))
# SST: same steps as theta without the 0-4 clip; script for DataLL_SST not published (LOG A9)
S = SS.astype(float); S[np.isnan(S).sum(1) > 224] = np.nan
k = np.isfinite(S).any(1); F = S[k]
for j in range(F.shape[0]): F[j] = lib.fill_linear(F[j])
S[k] = lib.lowess_uniform(F.T, 5).T
ok = np.isfinite(TH).all(1); okS = np.isfinite(S).all(1); both = ok & okS
print('n theta', ok.sum(), 'n sst', okS.sum(), 'both', both.sum())
R = {}
def rep(k, v, tgt=''):
    R[k] = [float(v), tgt]; print(f'{k:44s}{v:12.5g}   paper: {tgt}')
masks = {'nomask': ok}
if len(sys.argv) > 3:
    masks['mask'] = ok & np.load(sys.argv[3]).ravel('F')
out = {}
for nm, m in masks.items():
    print('=====', nm, m.sum())
    ssy, ssm, st, by, bm = lib.anova_year_month(TH[m], yr20, month)
    rep(nm+' grid points', m.sum(), '20822')
    rep(nm+' median seasonal SS frac', np.median(ssm/st), '0.59 ("59% of the variance")')
    rep(nm+' median interannual SS frac', np.median(ssy/st), '-')
    rep(nm+' area-wt mean seasonal frac', np.average(ssm/st, weights=W[m]), '-')
    D = lib.detrend(TH[m], timenum)
    maps, pcs, ev = lib.eof(D, 8)
    rep(nm+' EOF1-4 %', ev[:4].sum(), '59.4'); rep(nm+' EOF5-8 %', ev[4:].sum(), '4.2'); rep(nm+' EOF1-8 %', ev.sum(), '63.6')
    print('  each', ev.round(2))
    T = lib.trend(TH[m], timenum)*365.25
    rep(nm+' median theta trend /yr', np.median(T), '-1.3e-3')
    rep(nm+' area-wt frac theta trend<0', W[m][T < 0].sum()/W[m].sum(), '-')
    mS = m & okS; TS = lib.trend(S[mS], timenum)*365.25
    rep(nm+' SST warming area frac', W[mS][TS > 0].sum()/W[mS].sum(), '0.89')
    rep(nm+' SST warming frac, unweighted', (TS > 0).mean(), '0.89')
    rep(nm+' median SST trend C/yr', np.median(TS), '2.3e-2')
    rep(nm+' SST change over period, C', np.average(TS, weights=W[mS])*(time[-1]-time[0]), '0.43')
    _, _, evS = lib.eof(lib.detrend(S[mS], timenum), 8)
    rep(nm+' EOF SST 1-8 %', evS.sum(), '89.7'); rep(nm+' EOF SST ratio (5-8)/(1-4)', evS[4:].sum()/evS[:4].sum(), '~1:4 claimed')
    rep(nm+' EOF theta ratio (5-8)/(1-4)', ev[4:].sum()/ev[:4].sum(), '~1:4 claimed')
    # fig S16
    _, _, _, byS, bmS = lib.anova_year_month(S[mS], yr20, month)
    idx = np.isin(np.where(m)[0], np.where(mS)[0])
    rep(nm+' r season dSST vs dtheta', np.corrcoef(bmS.ravel(), bm[:, idx].ravel())[0, 1], '-0.51')
    rep(nm+' r annual dSST vs dtheta', np.corrcoef(byS.ravel(), by[:, idx].ravel())[0, 1], '-0.23')
    rep(nm+' N season', bmS.size, '354240'); rep(nm+' N annual', byS.size, '590400')
    rep(nm+' r trend theta vs trend SST', np.corrcoef(T[idx], TS)[0, 1], 'not in text')
    np.savez(f'{o}/sat_{nm}.npz', m=m, T=T, mS=mS, TS=TS, pcs=pcs, ev=ev, maps=maps, by=by, bm=bm, ssy=ssy, ssm=ssm, st=st, LAT=LAT, LON=LON, time=time, gm=np.exp(np.nanmean(np.log(np.where(TH[m] > 0, TH[m], np.nan)), 1)), mean=TH[m].mean(1))
json.dump(R, open(f'{o}/sat_results.json', 'w'), indent=1)
