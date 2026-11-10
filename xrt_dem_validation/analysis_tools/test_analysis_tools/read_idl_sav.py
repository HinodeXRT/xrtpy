from pathlib import Path

import numpy as np
from scipy.io import readsav

sav_path = Path("idl_dem_idx159_10Oct2008_7filters_MC0.sav")
data = readsav(str(sav_path), python_dict=True)

logT    = np.array(data['logt_out']).ravel()
dem_idl = np.array(data['dem_out']).ravel()   # 1D for MC0 case
chisq   = float(data['chisq'])
filters = [f.decode() if isinstance(f, bytes) else f for f in data['obs_index']]
obs_val = np.array(data['obs_val'])
obs_date = data['observation_date'].decode() if isinstance(data['observation_date'], bytes) else data['observation_date']

print(f"Date:       {obs_date}")
print(f"Filters:    {filters}")
print(f"Intensities:{obs_val}")
print(f"logT range: {logT.min():.2f} to {logT.max():.2f}  n={len(logT)}")
print(f"DEM shape:  {dem_idl.shape}")
print(f"DEM[:5]:    {dem_idl[:5]}")
print(f"chi2:       {chisq:.6f}")

print("\n--- Filters stored in SAV ---")
filters_sav = [f.decode() if isinstance(f, bytes) else f for f in data['obs_index']]
print(filters_sav)

print("\n--- Intensities stored in SAV ---")
print(np.array(data['obs_val']))
