"""Ocean mask substitute for the authors' unpublished BathyLL.mat.
Mean ETOPO5 depth per 1-degree cell < -300 m, on the 82x360 low-latitude grid (40.5N-40.5S, lon 20.5-379.5).
Input: global_topography.mat from the Climate Data Toolbox (github.com/chadagreene/CDT, cdt/cdt_data/).
usage: build_mask.py <global_topography.mat> <out.npy>"""
import sys, numpy as np, scipy.io as sio
m = sio.loadmat(sys.argv[1]); Z = m['Z'].astype(float); lat = m['lat'].ravel(); lon = m['lon'].ravel()
Z = Z[np.argsort(-lat)]                                  # north to south, 1/12 degree
mean = Z.reshape(180, 12, 360, 12).mean((1, 3))          # 1-degree cell means
first = np.floor(lon[0]); O2 = 20.5+np.arange(360.)
idx = ((O2-0.5-first) % 360).astype(int)
np.save(sys.argv[2], mean[49:131][:, idx] < -300)
