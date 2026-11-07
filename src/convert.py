"""Convert the Zenodo csv to npy. Layout found: rows 0-1 Year/Day, then 64,800 rows, 897 time columns.
Col 0/1 labels are lat/lon vectors recycled, so they are not a true lat-lon pair per row (LOG D12)."""
import sys, numpy as np, pandas as pd
src, dst = sys.argv[1:3]
h = pd.read_csv(src, header=None, nrows=2).values[:, 2:].astype(int)
b = pd.read_csv(src, header=None, skiprows=2, dtype=np.float32, engine='c').values
print(b.shape, np.unique(b[:,0]).size, np.unique(b[:,1]).size)
np.save(dst+'_time.npy', h); np.save(dst+'_labels.npy', b[:,:2]); np.save(dst+'.npy', b[:,2:])
