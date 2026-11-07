"""Python equivalents of the MATLAB steps in the authors' repo (ljustick/genomic_to_space_nut_stress v1.0.2).
Each function names the script/function it stands in for."""
import numpy as np

def daylength(yd, lat):
    """calc_thetaprime_v2.m daylength(). Hours of daylight."""
    dec = 23.5*np.cos(2*np.pi*(yd-172)/365)
    c = np.clip(-np.tan(np.radians(lat))*np.tan(np.radians(dec)), -1, 1)
    return np.degrees(np.arccos(c))/15*2

def bbp443_to_bbp532(bbp443, sbp):
    """Power-law conversion named in the code comment; the wavelength form is an assumption (see LOG A3)."""
    return bbp443*(443/532)**sbp

def thetaprime(chl, bbp532, k490, par, mld, yd, lat):
    """calc_thetaprime_v2.m. Fixes the repo typo `chlCSM` -> CChlSM (LOG D1)."""
    cphy = (bbp532-0.00029)*13000
    cphy = np.where(cphy < 0, np.nan, cphy)
    obs = cphy/chl
    kpar = 0.0665+0.874*k490-0.00121/k490
    pardl = par/daylength(yd, lat)
    ig = 0.98*pardl*np.exp(-kpar*(mld/2))
    sm = (1+np.exp(-0.15*pardl))/(1+np.exp(-3*ig))
    return 150*sm/obs

def fill_linear(y):
    """Stand-in for inpaint_nans(y,1) on a vector: linear in gaps, linear extrapolation at the ends (LOG A5)."""
    y = y.copy(); ok = np.isfinite(y)
    if ok.all() or ok.sum() < 2: return y
    x = np.arange(y.size); xo = x[ok]; yo = y[ok]
    out = np.interp(x, xo, yo)
    lo = x < xo[0]; hi = x > xo[-1]
    out[lo] = yo[0]+(x[lo]-xo[0])*(yo[1]-yo[0])/(xo[1]-xo[0])
    out[hi] = yo[-1]+(x[hi]-xo[-1])*(yo[-1]-yo[-2])/(xo[-1]-xo[-2])
    return out

def lowess_uniform(Y, span):
    """MATLAB smooth(y,span,'lowess') for evenly spaced points, along axis 0. Y: (nt, ...).
    Even spans are reduced to the next odd number (LOG A6)."""
    Y = np.asarray(Y, float); shape = Y.shape; Y = Y.reshape(shape[0], -1)
    n = Y.shape[0]; span = min(span, n); span -= (span % 2 == 0)
    h = (span-1)//2
    if h == 0: return Y.reshape(shape)
    d = np.abs(np.arange(-h, h+1)).astype(float); w = (1-(d/h)**3)**1.5
    # interior: symmetric window, local line reduces to weighted mean
    out = np.empty_like(Y)
    W = w/w.sum()
    for c in range(Y.shape[1]):
        out[h:n-h, c] = np.convolve(Y[:, c], W, mode='valid')
    # ends: weighted line through the first/last `span` points
    x = np.arange(span, dtype=float)
    for i in range(h):
        for idx, tgt, seg in ((i, x[i], np.s_[:span]), (n-1-i, x[span-1-i], np.s_[n-span:])):
            dd = np.abs(x-tgt); dm = dd.max(); ww = np.sqrt(np.clip(1-(dd/dm)**3, 0, None)**3)
            A = np.c_[np.ones(span), x-tgt]
            coef = np.linalg.lstsq(A*ww[:, None], Y[seg]*ww[:, None], rcond=None)[0]
            out[idx] = coef[0]
    return out.reshape(shape)

def preprocess(D, ocean_mask, nmax=4.0, nan_frac=0.25):
    """preprocess_LowLat_thetaPrime_v3_bathy.m. D: (npix, nt). >25% NaN dropped (C>224 of 897),
    clip [0,Nmax], gap fill, 5-point lowess, clip again."""
    D = D.astype(float).copy(); nt = D.shape[1]
    D[~ocean_mask] = np.nan
    limit = 224 if nt == 897 else int(np.floor(nan_frac*nt))
    D[np.isnan(D).sum(1) > limit] = np.nan
    D = np.clip(D, 0, nmax)
    keep = np.isfinite(D).any(1)
    F = D[keep]
    for j in range(F.shape[0]): F[j] = fill_linear(F[j])
    S = lowess_uniform(F.T, 5).T
    D[keep] = np.clip(S, 0, nmax)
    return D

def trend(D, t):
    """CDT trend(): OLS slope per unit t, for rows of D (npix, nt), complete rows only."""
    tc = t-t.mean()
    return (D-D.mean(1, keepdims=True))@tc/(tc@tc)

def detrend(D, t):
    tc = t-t.mean(); b = trend(D, t)
    return D-D.mean(1, keepdims=True)-np.outer(b, tc)

def eof(D, n=8):
    """CDT eof(): unweighted SVD of time-anomaly data (rows pixels). Returns maps, pcs, % variance."""
    X = D-D.mean(1, keepdims=True)
    C = X.T@X
    ev, U = np.linalg.eigh(C); o = np.argsort(ev)[::-1]; ev = ev[o]; U = U[:, o]
    s = np.sqrt(np.clip(ev, 0, None))
    maps = (X@U[:, :n])/s[:n]
    pcs = U[:, :n]*s[:n]
    return maps, pcs, 100*ev[:n]/ev.sum()

def anova_year_month(D, years, months):
    """ANOVA_analysis.m: anovan(y,{year,month}) main effects, Type III SS (default), sum-to-zero coefficients.
    Returns SS_year, SS_month, SS_total, year effects, month effects."""
    yu = np.unique(years); mu = np.unique(months)
    def code(v, u):
        Z = np.zeros((v.size, u.size-1))
        for i, x in enumerate(u[:-1]): Z[v == x, i] = 1
        Z[v == u[-1]] = -1
        return Z
    Zy, Zm = code(years, yu), code(months, mu); one = np.ones((years.size, 1))
    full = np.c_[one, Zy, Zm]
    def rss(M):
        Q, _ = np.linalg.qr(M); r = D.T-Q@(Q.T@D.T); return (r**2).sum(0)
    R = rss(full)
    ssy, ssm = rss(np.c_[one, Zm])-R, rss(np.c_[one, Zy])-R
    sst = ((D-D.mean(1, keepdims=True))**2).sum(1)
    beta = np.linalg.lstsq(full, D.T, rcond=None)[0]
    by = beta[1:yu.size]; by = np.r_[by, -by.sum(0, keepdims=True)]
    bm = beta[yu.size:]; bm = np.r_[bm, -bm.sum(0, keepdims=True)]
    return ssy, ssm, sst, by, bm

def boot_ci_mean(x, n=10000, seed=1):
    """Percentile CI. MATLAB bootci defaults to BCa, so this is a stated approximation (LOG A7)."""
    r = np.random.default_rng(seed); m = np.empty(n)
    for i in range(0, n, 500):
        k = min(500, n-i); m[i:i+k] = x[r.integers(0, x.size, (k, x.size))].mean(1)
    return x.mean(), *np.percentile(m, [2.5, 97.5])
