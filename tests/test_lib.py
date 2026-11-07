import sys, numpy as np; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__file__), '..', 'src'))
import lib
r = np.random.default_rng(0)
# line is preserved by lowess
y = 2+0.3*np.arange(50.); assert np.allclose(lib.lowess_uniform(y, 5), y)
# gaps
g = y.copy(); g[[0, 7, 8, 49]] = np.nan; assert np.allclose(lib.fill_linear(g), y)
# trend/detrend
t = np.arange(100.); D = np.vstack([3*t+1, -t]); assert np.allclose(lib.trend(D, t), [3, -1]); assert np.allclose(lib.detrend(D, t), 0)
# EOF recovers two planted modes and sums to variance
T = np.arange(240.); a, b = np.sin(2*np.pi*T/24), np.cos(2*np.pi*T/60)
m1, m2 = r.normal(size=500), r.normal(size=500)
m, p, v = lib.eof(np.outer(m1, a)+np.outer(m2, b), 2); assert abs(v.sum()-100) < 1e-6
assert abs(np.corrcoef(p[:, 0], a)[0, 1]) > .99 or abs(np.corrcoef(p[:, 1], a)[0, 1]) > .99
# ANOVA agrees with statsmodels Type III on an unbalanced design
import pandas as pd, statsmodels.api as sm; from statsmodels.formula.api import ols
yr = np.repeat(np.arange(2002, 2012), 12)[:-5]; mo = np.tile(np.arange(1, 13), 10)[:-5]
Y = r.normal(size=yr.size)+.3*(yr-2006)+np.sin(mo)
dfr = pd.DataFrame(dict(Y=Y, y=yr, m=mo))
f = ols('Y~C(y,Sum)+C(m,Sum)', dfr).fit(); tab = sm.stats.anova_lm(f, typ=3)
sy, sm_, st, by, bm = lib.anova_year_month(Y[None], yr, mo)
assert np.isclose(sy[0], tab.loc['C(y, Sum)', 'sum_sq']) and np.isclose(sm_[0], tab.loc['C(m, Sum)', 'sum_sq'])
assert abs(by.sum()) < 1e-9 and abs(bm.sum()) < 1e-9
# theta prime hand value: chl .1, bbp .002, k490 .03, par 40, mld 50, yd 172, lat 0
tp = lib.thetaprime(.1, .002, .03, 40., 50., 172, 0.)
dl = 12.; k = .0665+.874*.03-.00121/.03; pdl = 40/dl; ig = .98*pdl*np.exp(-k*25)
assert np.isclose(tp, 150*(1+np.exp(-.15*pdl))/(1+np.exp(-3*ig))/((.002-.00029)*13000/.1))
print('all tests passed')
