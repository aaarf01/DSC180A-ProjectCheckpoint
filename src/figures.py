"""Rebuild the paper's main figures from the replication output. Static PNGs on a light surface."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.colors as mc, lib, sys
from scipy import stats
d, o = '../data', '../out'; F = f'{o}/fig'; import os; os.makedirs(F, exist_ok=True)
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#888', 'axes.facecolor': '#fcfcfb', 'figure.facecolor': '#fcfcfb', 'font.family': 'DejaVu Sans'})
BL, OR, GR = '#2a78d6', '#eb6834', '#8a8a86'
DIV = mc.LinearSegmentedColormap.from_list('d', ['#1c4f94', '#2a78d6', '#a9c8ee', '#f0efec', '#f3b3a6', '#d6432f', '#8c1d12'])
SEQ = mc.LinearSegmentedColormap.from_list('s', ['#e8f1fb', '#a9c8ee', '#2a78d6', '#123a70'])
z = np.load(f'{o}/sat_mask.npz'); m = z['m']; mS = z['mS']; time = z['time']
g = lambda v, mask=m: (lambda a: (a.__setitem__(mask, v), a.reshape(82, 360, order='F'))[1])(np.full(29520, np.nan))
EXT = [20, 380, -41, 41]
def mp(ax, A, vmin, vmax, cmap, title=''):
    A = np.ma.masked_invalid(A); c = plt.cm.get_cmap(cmap) if isinstance(cmap, str) else cmap
    c = c.copy(); c.set_bad('#d9d8d3')
    im = ax.imshow(A, extent=EXT, vmin=vmin, vmax=vmax, cmap=c, aspect='auto', interpolation='nearest')
    ax.set_xticks([45, 135, 225, 315]); ax.set_xticklabels(['45E', '135E', '135W', '45W']); ax.set_yticks([-40, 0, 40]); ax.set_title(title, loc='left', fontsize=9); return im
# Fig 1
x = pd.read_csv(f'{d}/data_insitu_v2_25-04.csv'); t = x.theta_prime_2x2degree_40days_avg; k = t.notna()
fig, ax = plt.subplots(1, 3, figsize=(10, 3.1), gridspec_kw=dict(width_ratios=[1.3, .7, 1.1]))
ax[0].scatter(x.insitu_and_woa_nutricline[k], t[k], s=6, c=BL, alpha=.35, lw=0)
b = pd.cut(x.insitu_and_woa_nutricline[k], np.arange(0, 241, 25)); md = t[k].groupby(b, observed=True).median()
ax[0].plot([i.mid for i in md.index], md.values, color=OR, lw=2, marker='o', ms=4, label='Bin median'); ax[0].set_ylim(0, 4.5); ax[0].legend(frameon=False)
ax[0].set(xlabel='Nutricline depth (m)', ylabel="Θ′", title='A  Stations only (paper uses all grid boxes)')
r = pd.read_csv(f'{o}/fig1B_r2.csv'); bars = ax[1].bar(['ZNO3 +\nbiomarkers', 'ZNO3'], r.r2.values[::-1], color=BL, width=.55)
for q, v in zip(bars, r.r2.values[::-1]): ax[1].text(q.get_x()+.27, v+.01, f'{v:.2f}', ha='center')
ax[1].set(ylabel='Out-of-bag R²', ylim=(0, .75), title='B  Random forest')
p = np.load(f'{o}/figS4_published.npy'); f = pd.read_csv(f'{o}/figS4_r2_asrepo.csv', index_col=0).values
ax[2].plot([-10, 70], [-10, 70], color=GR, lw=1); ax[2].scatter(p.ravel(), f.ravel(), s=10, c=BL); ax[2].set(xlabel='Published Fig. S4 R² (%)', ylabel='Replicated R² (%)', title='C  All 90 Fig. S4 cells')
plt.tight_layout(); plt.savefig(f'{F}/fig1.png', dpi=130); plt.close()
# Fig S4 true 1P
ff = pd.read_csv(f'{o}/figS4_r2_fixed.csv', index_col=0).values; D = np.arange(2, 31, 2)
fig, a = plt.subplots(figsize=(5, 2.8)); a.plot(D, p[:, 0], color=OR, lw=2, marker='o', ms=4, label='Published "1P" column'); a.plot(D, ff[:, 3], color=GR, lw=2, label='1° file (replicated)'); a.plot(D, ff[:, 0], color=BL, lw=2, marker='o', ms=4, label='Actual 1-pixel file (replicated)')
a.set(xlabel='Time window (± days)', ylabel='R² (%)'); a.legend(frameon=False); plt.tight_layout(); plt.savefig(f'{F}/figS4_bug.png', dpi=130); plt.close()
# Fig 2
bm = z['bm']; fig, ax = plt.subplots(3, 2, figsize=(10, 6.3))
im = mp(ax[0, 0], g(z['gm']), 0, 4, SEQ, 'A  Geometric mean Θ′'); plt.colorbar(im, ax=ax[0, 0], fraction=.03)
ax[0, 1].axis('off'); ax[0, 1].text(0, .5, 'B  WOA 2023 nutricline map:\nnot reproduced (NOAA server down)', color='#555')
for a, i, n in zip(ax.ravel()[2:], (0, 3, 6, 9), ('C  January', 'D  April', 'E  July', 'F  October')):
    im = mp(a, g(bm[i]), -1, 1, DIV, n+' anomaly'); plt.colorbar(im, ax=a, fraction=.03)
plt.tight_layout(); plt.savefig(f'{F}/fig2.png', dpi=130); plt.close()
# Fig 3
oni = np.load(f'{d}/oni.npy'); pdo = np.load(f'{d}/pdo.npy'); pdo = pdo[pdo[:, 1] < 90]
Y = lib.lowess_uniform(z['pcs'][:, 4:6], 30); by = z['by']
xo = np.interp(time, oni[:, 0], oni[:, 1]); xp = np.interp(time, pdo[:, 0], pdo[:, 1])
sg = np.sign(np.corrcoef(Y[:, 0], xo)[0, 1]), np.sign(np.corrcoef(Y[:, 1], xp)[0, 1])
fig, ax = plt.subplots(2, 3, figsize=(11, 4.6), gridspec_kw=dict(width_ratios=[1.25, 1, 1]))
for a, y, xx, nm, s in ((ax[0, 0], Y[:, 0]*sg[0], xo, 'ONI', 'A  PC5'), (ax[1, 0], Y[:, 1]*sg[1], xp, 'PDO', 'D  PC6')):
    zs = lambda v: (v-v.mean())/v.std()
    a.plot(time, zs(y), color=BL, lw=1.6, label=s[3:]); a.plot(time, zs(xx), color=OR, lw=1.6, label=nm); a.legend(frameon=False, ncol=2, loc='upper left'); a.set_xticks(range(2002,2023,4)); a.set_title(f'{s} vs {nm}, r = {abs(np.corrcoef(y, xx)[0, 1]):.2f} (both z-scored)', loc='left', fontsize=9)
yrs = np.arange(2002, 2022); Yr = lambda L: by[np.isin(yrs, L)].mean(0)
mp(ax[0, 1], g(by[yrs == 2011][0]), -1, 1, DIV, 'B  2011 (La Niña)'); mp(ax[0, 2], g(by[yrs == 2015][0]), -1, 1, DIV, 'C  2015 (El Niño)')
mp(ax[1, 1], g(Yr(list(range(2007, 2015))+[2020, 2021])), -.5, .5, DIV, 'E  PDO cool years'); mp(ax[1, 2], g(Yr(list(range(2002, 2007))+list(range(2015, 2020)))), -.5, .5, DIV, 'F  PDO warm years')
plt.tight_layout(); plt.savefig(f'{F}/fig3.png', dpi=130); plt.close()
# Fig 4
T, TS = z['T'], z['TS']
fig, ax = plt.subplots(1, 3, figsize=(11, 2.9), gridspec_kw=dict(width_ratios=[1.2, 1.2, .9]))
im = mp(ax[0], g(TS, mS), -.1, .1, DIV, 'A  SST trend (°C/yr)'); plt.colorbar(im, ax=ax[0], fraction=.03)
im = mp(ax[1], g(T), -.05, .05, DIV, 'B  Θ′ trend (per yr)'); plt.colorbar(im, ax=ax[1], fraction=.03)
ax[2].hist(T, bins=np.linspace(-.1, .1, 80), color=OR, alpha=.8, label='Θ′'); ax[2].hist(TS, bins=np.linspace(-.1, .1, 80), color=BL, alpha=.7, label='SST'); ax[2].axvline(0, color='#555', lw=1); ax[2].legend(frameon=False); ax[2].set(title='C  Trend histograms', xlabel='per year', yticks=[])
plt.tight_layout(); plt.savefig(f'{F}/fig4.png', dpi=130); plt.close()
# S10
fig, ax = plt.subplots(1, 3, figsize=(11, 2.8)); st = z['st']
im = mp(ax[0], g(z['ssy']/st), 0, .2, SEQ, 'A  Interannual share of variance'); plt.colorbar(im, ax=ax[0], fraction=.03)
im = mp(ax[1], g(z['ssm']/st), 0, 1, SEQ, 'B  Seasonal share of variance'); plt.colorbar(im, ax=ax[1], fraction=.03)
ax[2].hist(z['ssy']/st, bins=50, range=(0, 1), color=OR, label='Interannual'); ax[2].hist(z['ssm']/st, bins=50, range=(0, 1), color=BL, alpha=.75, label='Seasonal'); ax[2].legend(frameon=False); ax[2].set(title='C  Across grid points', xlabel='Share of variance')
plt.tight_layout(); plt.savefig(f'{F}/figS10.png', dpi=130); plt.close()
print('r PC5/ONI', np.corrcoef(Y[:, 0], xo)[0, 1], 'r PC6/PDO', np.corrcoef(Y[:, 1], xp)[0, 1])
