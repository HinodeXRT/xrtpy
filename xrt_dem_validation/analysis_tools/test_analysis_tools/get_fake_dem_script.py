# read_idl_sav.py
# Load and inspect an IDL .sav file from fake_dem_base_only.pro

from pathlib import Path

import numpy as np
from scipy.io import readsav

sav_path = Path("Fake_IDL_data/idl_dem_idx159_10Oct2008_7filters_MC0.sav")

data = readsav(str(sav_path), python_dict=True)

print("Keys in .sav file:")
for k, v in data.items():
    if hasattr(v, 'shape'):
        print(f"  {k:25s}  shape={v.shape}  dtype={v.dtype}")
    else:
        print(f"  {k:25s}  value={v}")

print("\n--- logT_out ---")
logT = np.array(data['logt_out']).ravel()
print(f"  shape: {logT.shape}")
print(f"  range: {logT.min():.2f} to {logT.max():.2f}")

print("\n--- dem_out ---")
dem = np.array(data['dem_out'])
print(f"  shape: {dem.shape}")
print(f"  dem_out[*,0] (base DEM) first 5: {dem[:,0][:5]}")

print("\n--- obs_index (filters) ---")
print(f"  {data.get('obs_index', 'NOT FOUND')}")

print("\n--- obs_val (intensities) ---")
print(f"  {data.get('obs_val', 'NOT FOUND')}")

print("\n--- observation_date ---")
print(f"  {data.get('observation_date', 'NOT FOUND')}")

print("\n--- chisq ---")
print(f"  {data.get('chisq', 'NOT FOUND')}")
