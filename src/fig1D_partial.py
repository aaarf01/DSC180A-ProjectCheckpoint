"""Partial Fig. 1D: correlation of each N-trait gene with satellite theta' (2 deg, +-20 d) across samples.
Gene frequencies come from Ustick et al. 2021 data S1 (Science supplement, not redistributed here; download it yourself).
Trait-to-gene mapping from Martiny et al. 2026 table S1. Genes missing from data S1 are reported as unavailable.
usage: fig1D_partial.py <abe6301-ustick-sm-data-s1.csv> <data_insitu_v2_25-04.csv> <out dir>"""
import sys, numpy as np, pandas as pd
from scipy import stats
s1p, insitu, o = sys.argv[1:4]
s1 = pd.read_csv(s1p); z = pd.read_csv(insitu)
s1['key'] = s1.Sample.str.upper(); z['key'] = z.Sample.str.upper()
m = s1.merge(z[['key', 'theta_prime_2x2degree_40days_avg']], on='key')
t = m.theta_prime_2x2degree_40days_avg; k = t.notna()
traits = {'Nitrile': ['nitC', 'nitD', 'hypA', 'hypB'], 'Cyanate': ['cynA', 'cynB', 'cynD', 'cynS'],
          'Nitrate': ['narB', 'napA', 'nirA2', 'narX'], 'Urea': ['ureA', 'ureB', 'ureC', 'ureD', 'ureE', 'ureF', 'urtA', 'speB'],
          'Nitrite': ['nirA', 'focA'], 'Ammonia': ['amt1', 'amt']}
rows = []
for tr, gs in traits.items():
    for g in gs:
        if g in m.columns:
            f = m.loc[k, g]; r, p = stats.pearsonr((f-f.mean())/f.std(), t[k]); rows.append((tr, g, r, p, int(k.sum()), 'available'))
        else:
            rows.append((tr, g, np.nan, np.nan, 0, 'not in published gene table'))
R = pd.DataFrame(rows, columns=['trait', 'gene', 'r', 'p', 'n', 'status']); R.to_csv(f'{o}/fig1D_partial.csv', index=False)
print(R.round(3).to_string(index=False))
