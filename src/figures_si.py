"""Supplementary figures that the public data supports: S5, S7, S8/S9, S11, S12, S13, S14, S15, S16, S17."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, lib
from scipy.signal import periodogram
from scipy import stats
exec(open('figures.py').read().split('# Fig 1\n')[0])
S=None
by, bm, pcs, ev, maps = z['by'], z['bm'], z['pcs'], z['ev'], z['maps']
mon = 'Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec'.split()
fig, ax = plt.subplots(4, 3, figsize=(10, 6.6))
for i, a in enumerate(ax.ravel()): mp(a, g(bm[i]), -1, 1, DIV, mon[i]); a.set_xticks([]); a.set_yticks([])
plt.tight_layout(); plt.savefig(f'{F}/figS7.png', dpi=110); plt.close()
fig, ax = plt.subplots(8, 3, figsize=(11, 10), gridspec_kw=dict(width_ratios=[2.2, 1, 1.2]))
per = []
for i in range(8):
    ax[i, 0].plot(time, pcs[:, i], color=BL, lw=.8); ax[i, 0].set_ylabel(f'PC{i+1}\n{ev[i]:.1f}%')
    f, P = periodogram(pcs[:, i]-pcs[:, i].mean(), fs=365.25/8); ax[i, 1].plot(f, P/P.max(), color=BL, lw=1); ax[i, 1].set_xlim(0, 3); ax[i, 1].set_yticks([])
    per.append(f[1:][np.argmax(P[1:])]); mp(ax[i, 2], g(maps[:, i]), -.05, .05, DIV, ''); ax[i, 2].set_xticks([]); ax[i, 2].set_yticks([])
ax[0, 1].set_title('Power (cycles/yr)', fontsize=9); ax[0, 0].set_title('Temporal mode', fontsize=9); ax[0, 2].set_title('Spatial loading', fontsize=9)
plt.tight_layout(); plt.savefig(f'{F}/figS8.png', dpi=100); plt.close()
print('PC peak cycles/yr', np.round(per, 2))
fig, ax = plt.subplots(4, 5, figsize=(11, 4.6)); yrs = range(2002, 2022)
for i, a in enumerate(ax.ravel()): mp(a, g(by[i]), -1, 1, DIV, str(list(yrs)[i])); a.set_xticks([]); a.set_yticks([])
plt.tight_layout(); plt.savefig(f'{F}/figS14.png', dpi=110); plt.close()
TH = np.load(f'{d}/theta_LL.npy').astype(float)
t2 = np.load(f'{d}/time.npz'); mo = t2['month']
fig, ax = plt.subplots(1, 4, figsize=(11, 1.9))
for a, ms, n in zip(ax, ((12, 1, 2), (3, 4, 5), (6, 7, 8), (9, 10, 11)), ('Winter DJF', 'Spring MAM', 'Summer JJA', 'Fall SON')):
    mp(a, g(np.nanmean(TH[m][:, np.isin(mo, ms)], 1)), 0, 4, SEQ, n); a.set_xticks([]); a.set_yticks([])
plt.tight_layout(); plt.savefig(f'{F}/figS11.png', dpi=110); plt.close()
LAT, LON = z['LAT'], z['LON']
fig, ax = plt.subplots(5, 1, figsize=(8, 6), sharex=True)
for a, (n, la_, lo_) in zip(ax, (('HOT', 23, 202), ('BATS', 32, 296), ('SPac', -25, 230), ('EqAtl', 0, 350), ('EqPac', 1, 170))):
    i = np.argmin(np.hypot(LAT-la_, (LON-lo_+180) % 360-180)+1e3*(~m)); y = TH[i]; rec = (np.where(m)[0] == i)
    r = y.mean()+(maps[rec][0]*pcs).sum(1); a.plot(time, y, color=GR, lw=.8, label='Θ′'); a.plot(time, r, color=BL, lw=1, label='Top 8 modes'); a.set_ylabel(n)
ax[0].legend(frameon=False, ncol=2); plt.tight_layout(); plt.savefig(f'{F}/figS13.png', dpi=110); plt.close()
# S5, S17 from station data
r = pd.read_csv(f'{o}/rf_oob.csv'); fit = stats.linregress(r.nitri, r.theta); an = r.theta-(fit.intercept+fit.slope*r.nitri)
print('linear fit', fit.slope, fit.intercept)
print('anomaly vs omega', {c: round(stats.pearsonr(an, r[c])[0], 3) for c in ['p_lim_h', 'n_lim_h', 'n_lim_m', 'fe_lim_h', 'fe_lim_m']})
fig, ax = plt.subplots(1, 2, figsize=(9, 2.7)); b = np.linspace(-3, 3, 50)
ax[0].hist(r.theta-r.nut_pred, bins=b, color=OR, alpha=.8, label='Nutricline only'); ax[0].hist(r.theta-r.full_pred, bins=b, color=BL, alpha=.7, label='Nutricline + biomarkers'); ax[0].legend(frameon=False); ax[0].set(xlabel='Observed minus predicted Θ′', yticks=[], title='S17A  Residuals')
ax[1].scatter(r.theta-r.nut_pred, r.theta-r.full_pred, c=r.p_lim_h, s=7, cmap=SEQ); ax[1].set(xlabel='Residual, nutricline only', ylabel='Residual, full', xlim=(-3, 4), ylim=(-3, 4), title='S17B  Colored by high P biomarker')
plt.tight_layout(); plt.savefig(f'{F}/figS17.png', dpi=110); plt.close()
