from pathlib import Path
from utils_case_io import read_mc_intensities_csv,run_dem_for_mc_csv, load_idl_dem_sav, case_dir


case_dir = case_dir("case_20080104_110426")
csv_path = case_dir / "mc_intensities_20080104_110426_IDL.csv"
idl = load_idl_dem_sav(case_dir / "xrt_dem_output_20080104_110426_MCITER100.sav")
mc = read_mc_intensities_csv(csv_path)


print(mc.filters)
print(mc.mc_intensities.shape)  # (N, n_filters)
print(mc.df.head())

observation_date = "2008-01-04T11:04:26"

out = run_dem_for_mc_csv(
    csv_path=csv_path,
    observation_date=observation_date,
    # intensity_errors=None,  # keep None to match IDL default behavior 
)

#import pdb; pdb.set_trace()

print("filters:", out.filters)
print("logT shape:", out.logT.shape)
print("dem_runs shape:", out.dem_runs.shape)
print("modeled_runs shape:", out.modeled_runs.shape)
print("chisq_runs shape:", out.chisq_runs.shape)

    
    
print('\n\n-New Section-\n')
print("IDL logT:", idl.logT.shape)
print("IDL dem_base:", idl.dem_base.shape)
print("IDL dem_runs:", idl.dem_runs.shape, "n_runs=", idl.n_runs)

