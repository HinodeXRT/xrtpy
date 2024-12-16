from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
from astropy import units as u
#from astropy.utils.data import get_pkg_data_filenames

from xrtpy.response.channel import Channel
from xrtpy.response.effective_area import EffectiveAreaFundamental

channel_names = [
    "Al-mesh",
    "Al-poly",
    "C-poly",
    "Ti-poly",
    "Be-thin",
    "Be-med",
    "Al-med",
    "Al-thick",
    "Be-thick",
    "Al-poly/Al-mesh",
    "Al-poly/Ti-poly",
    "Al-poly/Al-thick",
    "Al-poly/Be-thick",
    "C-poly/Ti-poly",
]

channel_single_filter_names = [
    "Al-mesh",
    "Al-poly",
    "C-poly",
    "Ti-poly",
    "Be-thin",
    "Be-med",
    "Al-med",
    "Al-thick",
    "Be-thick",
]

valid_dates = [
    datetime(year=2006, month=9, day=25, hour=22, minute=1, second=1),
    datetime(year=2007, month=9, day=22, hour=22, minute=1, second=1),
    datetime(year=2009, month=9, day=22, hour=22, minute=1, second=1),
    datetime(year=2010, month=9, day=22, hour=22, minute=1, second=1),
    datetime(year=2012, month=9, day=22, hour=22, minute=1, second=1),
    datetime(year=2015, month=9, day=22, hour=22, minute=1, second=1),
    datetime(year=2017, month=9, day=22, hour=22, minute=1, second=1),
    datetime(year=2019, month=9, day=23, hour=22, minute=1, second=1),
    datetime(year=2020, month=9, day=23, hour=22, minute=1, second=1),
    datetime(year=2021, month=9, day=23, hour=22, minute=1, second=1),
    datetime(year=2022, month=9, day=23, hour=22, minute=1, second=1),
]

invalid_dates = [
    datetime(year=2006, month=8, day=25, hour=22, minute=1, second=1),
    datetime(year=2005, month=9, day=22, hour=22, minute=1, second=1),
    datetime(year=2002, month=9, day=22, hour=22, minute=1, second=1),
    datetime(year=2000, month=9, day=22, hour=22, minute=1, second=1),
    datetime(year=1990, month=9, day=22, hour=22, minute=1, second=1),
]


@pytest.mark.parametrize("channel_name", channel_names)
def test_channel_name(channel_name):
    channel = Channel(channel_name)
    assert channel.name == channel_name


@pytest.mark.parametrize("name", channel_names)
def test_EffectiveArea_filter_name(name):
    instance = EffectiveAreaFundamental(
        name, datetime(year=2013, month=9, day=22, hour=22, minute=0, second=0)
    )
    assert instance.name == name


@pytest.mark.parametrize("date", valid_dates)
@pytest.mark.parametrize("name", channel_names)
def test_EffectiveArea_contamination_on_CCD(name, date):
    instance = EffectiveAreaFundamental(name, date)
    assert 0 <= instance.contamination_on_CCD <= 1206


@pytest.mark.parametrize("date", valid_dates)
@pytest.mark.parametrize("name", channel_single_filter_names)
def test_EffectiveArea_contamination_on_filter(name, date):
    instance = EffectiveAreaFundamental(name, date)
    assert 0 <= instance.contamination_on_filter <= 2901


@pytest.mark.parametrize("date", invalid_dates)
@pytest.mark.parametrize("name", channel_names)
def test_EffectiveArea_exception_is_raised(name, date):
    with pytest.raises(ValueError, match="Invalid date"):
        EffectiveAreaFundamental(name, date)


# def get_IDL_data_files(): 
#     filter_data_files = []
#     for dir in get_pkg_data_filenames(
#         "data/effective_area_IDL_testing_files", package="xrtpy.response.tests" 
#     ):
#         print(dir)
#         filter_data_files += list(Path(dir).glob("*.txt"))
#     return sorted(filter_data_files)

def get_IDL_data_files():
    data_dir = Path(__file__).parent / "data" / "effective_area_IDL_testing_files" 
    assert data_dir.exists(), f"Data directory {data_dir} does not exist."
    files = sorted(data_dir.glob("**/*.txt")) 
    #print(f"\n\n\nFound files: {files}\n\n")  # Debugging output
    return files


# #Working testing
# @pytest.mark.parametrize("filename", get_IDL_data_files())
# def test_effective_area_compare_idl(filename):
#     print(f"\n\nTesting file: {filename}\n")
    
#     # Read the filter name and observation date from the file
#     with filename.open() as f:
#         filter_name = f.readline().split()[1]
#         filter_obs_date = " ".join(f.readline().split()[1:])
#         print(f"Filter name: {filter_name}, Observation date: {filter_obs_date}")  # Debugging output

#     # Correct non-standard date format
#     filter_obs_date = filter_obs_date.replace("Sept", "Sep")
    
#     # Load IDL data from the file
#     IDL_data = np.loadtxt(filename, skiprows=3)
#     IDL_wavelength = IDL_data[:, 0] * u.AA
#     IDL_effective_area = IDL_data[:, 1] * u.cm**2
#     #print(f"Loaded IDL data: {IDL_data.shape}\n")
    
#     # Compute effective area using XRTpy
#     instance = EffectiveAreaFundamental(filter_name, filter_obs_date)
#     actual_effective_area = instance.effective_area()
#     #print(f"Calculated effective area (first 5): {actual_effective_area[:5]}\n")
    
#     # # Interpolate IDL effective area values to align with XRTpy wavelengths
#     # IDL_effective_area_interp = np.interp(instance.wavelength, IDL_wavelength, IDL_effective_area)
#     # #print(f"Interpolated IDL effective area (first 5): {IDL_effective_area_interp[:5]}\n")
#         # # Compare XRTpy and IDL effective areas using rtol=1e-4
#     # assert u.allclose(
#     #     actual_effective_area,
#     #     IDL_effective_area_interp,
#     #     rtol=1e-1,
#     # ), f"Effective areas differ for filter {filter_name} on {filter_obs_date}"

#         XRTpy_effective_area = np.interp(
#             IDL_wavelength,           # Target grid (IDL wavelengths)
#             instance.wavelength,      # Source grid (XRTpy wavelengths)
#             actual_effective_area     # Data to interpolate
#         )

#         assert u.allclose(
#         XRTpy_effective_area,    # Interpolated XRTpy values
#         IDL_effective_area,      # Original IDL values
#         rtol=1e-4,               # Relative tolerance
#     ), f"Effective areas differ for filter {filter_name} on {filter_obs_date}"


################################################################################################
################################################################################################
################################################################################################
# ################################################################################################
################################################################################################
################################################################################################
# ################################################################################################
import numpy as np
import pytest
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
from astropy import units as u
from xrtpy.response.effective_area import EffectiveAreaFundamental

# Directory to save plots
OUTPUT_DIR = Path("effective_area_plots")
OUTPUT_DIR.mkdir(exist_ok=True)

# Store test data for each filter
filter_test_data = {}


@pytest.mark.parametrize("filename", get_IDL_data_files())
def test_effective_area_compare_idl(filename):
    """
    Compare effective area calculations between XRTpy and IDL.
    """
    print(f"\n\nTesting file: {filename}\n")
    
    # Read the filter name and observation date from the file
    with filename.open() as f:
        filter_name = f.readline().split()[1]
        filter_obs_date = " ".join(f.readline().split()[1:])
        print(f"Filter name: {filter_name}, Observation date: {filter_obs_date}")
    
    # Correct non-standard date format
    filter_obs_date = filter_obs_date.replace("Sept", "Sep")
    
    # Load IDL data
    IDL_data = np.loadtxt(filename, skiprows=3)
    IDL_wavelength = IDL_data[:, 0] * u.AA
    IDL_effective_area = IDL_data[:, 1] * u.cm**2
    
    # Compute effective area using XRTpy
    instance = EffectiveAreaFundamental(filter_name, filter_obs_date)
    actual_effective_area = instance.effective_area()
    
    # Interpolate XRTpy effective area onto the IDL wavelength grid
    XRTpy_effective_area = np.interp(
        IDL_wavelength.value,           # Target grid (IDL wavelengths)
        instance.wavelength.value,      # Source grid (XRTpy wavelengths)
        actual_effective_area.value     # Data to interpolate
    )
    
    # Compare using relative tolerance
    rtol = 1e-4
    differences = np.abs(XRTpy_effective_area - IDL_effective_area.value)
    failed_indices = np.where(differences > rtol * np.abs(IDL_effective_area.value))[0]
    
    # Store test data for combined plots
    if filter_name not in filter_test_data:
        filter_test_data[filter_name] = []
    filter_test_data[filter_name].append(
        (IDL_wavelength.value, IDL_effective_area.value, XRTpy_effective_area, failed_indices, filter_obs_date)
    )
    
    # Plot results if there are failures
    if failed_indices.size > 0:
        print(f"Test failed for filter {filter_name} on {filter_obs_date}")
        plot_effective_area_comparison(
            IDL_wavelength.value,
            IDL_effective_area.value,
            XRTpy_effective_area,
            failed_indices,
            filter_name,
            filter_obs_date,
        )
        assert False, f"Effective areas differ for filter {filter_name} on {filter_obs_date}"

def plot_effective_area_comparison(wavelength, IDL_area, XRTpy_area, failed_indices, filter_name, obs_date):
    """
    Plot individual effective area comparison with failed points highlighted.
    """
    filter_dir = OUTPUT_DIR / filter_name
    filter_dir.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(10, 6))
    plt.plot(wavelength, IDL_area, label="IDL Effective Area", color="blue", lw=2)
    plt.plot(wavelength, XRTpy_area, label="XRTpy Effective Area", color="green", linestyle="--", lw=2)
    plt.scatter(
        wavelength[failed_indices],
        IDL_area[failed_indices],
        color="red",
        label=f"Failed Points: {len(failed_indices)}",
        zorder=5,
    )
    plt.xlabel("Wavelength (Å)", fontsize=14)
    plt.ylabel("Effective Area (cm²)", fontsize=14)
    plt.title(
        f"Effective Area Comparison for {filter_name} on {obs_date}\n"
        f"IDL Linear Interpolation and XRTpy Wavelength Grid",
        fontsize=16
    )
    plt.legend(fontsize=12)
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.tight_layout()
    
    output_file = filter_dir / f"{filter_name}_{obs_date.replace(':', '-')}_XRTpy_Wavelength_Grid.png"
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Saved plot to {output_file}")


def save_combined_filter_plots(filter_name, test_data):
    """
    Save a combined plot with subplots for all test dates for a single filter.
    """
    # Create a subdirectory for the filter
    filter_dir = OUTPUT_DIR / filter_name
    filter_dir.mkdir(parents=True, exist_ok=True)
    
    # Define the output file
    output_file = filter_dir / f"{filter_name}_IDL_Linear_combined_XRTpy_Wavelength_Grid.png"
    
    # Determine subplot grid size
    num_tests = len(test_data)
    num_cols = 3
    num_rows = (num_tests + num_cols - 1) // num_cols  # Ceiling division
    
    # Create the figure and axes
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(15, 5 * num_rows), squeeze=False)
    fig.suptitle(
        f"Combined Effective Area Comparison for {filter_name}\n"
        f"IDL Linear Interpolation and XRTpy Wavelength Grid",
        fontsize=18,
        y=1.02  # Adjust vertical spacing for the title
    )
    
    for i, (wavelength, IDL_area, XRTpy_area, failed_indices, obs_date) in enumerate(test_data):
        row, col = divmod(i, num_cols)
        ax = axes[row, col]
        
        # Plot data
        ax.plot(wavelength, IDL_area, label="IDL Effective Area", color="blue", lw=2)
        ax.plot(wavelength, XRTpy_area, label="XRTpy Effective Area", color="green", linestyle="--", lw=2)
        ax.scatter(
            wavelength[failed_indices],
            IDL_area[failed_indices],
            color="red",
            label=f"Failed Points: {len(failed_indices)}",
            zorder=5,
        )
        
        # Customize subplot
        ax.set_title(f"{obs_date}", fontsize=10)
        ax.set_xlabel("Wavelength (Å)", fontsize=8)
        ax.set_ylabel("Effective Area (cm²)", fontsize=8)
        ax.grid(True, linestyle="--", alpha=0.7)
        ax.legend(fontsize=6)
    
    # Remove any empty subplots
    for j in range(i + 1, num_rows * num_cols):
        fig.delaxes(axes[j // num_cols, j % num_cols])
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved combined plot to {output_file}")
# Generate combined plots for all filters after tests
@pytest.fixture(scope="session", autouse=True)
def generate_combined_plots(request):
    """
    Generate combined plots for each filter after all tests.
    """
    @request.addfinalizer
    def save_plots():
        for filter_name, test_data in filter_test_data.items():
            save_combined_filter_plots(filter_name, test_data)

################################################################################################
################################################################################################
# ################################################################################################
################################################################################################
################################################################################################
# ################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
# import matplotlib.pyplot as plt
# from pathlib import Path
# from datetime import datetime
# import numpy as np
# import pytest
# from astropy import units as u

# # Directory to save plots
# OUTPUT_DIR = Path("test_plots")
# OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# # Dictionary to store data for each filter
# filter_test_data = {}

# @pytest.mark.parametrize("filename", get_IDL_data_files())
# def test_effective_area_compare_idl(filename):
#     # Read the filter name and observation date from the file
#     with filename.open() as f:
#         filter_name = f.readline().split()[1]
#         filter_obs_date = " ".join(f.readline().split()[1:])
#         filter_obs_date = filter_obs_date.replace("Sept", "Sep")
    
#     # Load IDL data from the file
#     IDL_data = np.loadtxt(filename, skiprows=3)
#     IDL_wavelength = IDL_data[:, 0] * u.AA
#     IDL_effective_area = IDL_data[:, 1] * u.cm**2

#     # Compute effective area using XRTpy
#     instance = EffectiveAreaFundamental(filter_name, filter_obs_date)
#     actual_effective_area = instance.effective_area()

#     # Interpolate XRTpy effective area onto the IDL wavelength grid
#     XRTpy_effective_area_interp = np.interp(
#         IDL_wavelength.value,
#         instance.wavelength.value,
#         actual_effective_area.value,
#     )

#     # Make both arrays unitless
#     IDL_effective_area_unitless = IDL_effective_area.value
#     XRTpy_effective_area_unitless = XRTpy_effective_area_interp

#     # Compare effective areas
#     rtol = 1e-4
#     differences = np.abs(XRTpy_effective_area_unitless - IDL_effective_area_unitless)
#     failed_indices = np.where(differences > rtol * np.abs(IDL_effective_area_unitless))[0]

#     # Store test data
#     if filter_name not in filter_test_data:
#         filter_test_data[filter_name] = []
#     filter_test_data[filter_name].append((
#         IDL_wavelength.value,
#         IDL_effective_area_unitless,
#         XRTpy_effective_area_unitless,
#         failed_indices,
#         filter_obs_date,
#     ))

#     # Plot individual tests
#     save_effective_area_plot(
#         IDL_wavelength.value,
#         IDL_effective_area_unitless,
#         XRTpy_effective_area_unitless,
#         failed_indices,
#         filter_name,
#         filter_obs_date,
#     )

#     # Fail test if there are mismatched points
#     if failed_indices.size > 0:
#         assert False, f"Effective areas differ for filter {filter_name} on {filter_obs_date}"


# def save_effective_area_plot(wavelength, IDL_area, XRTpy_area, failed_indices, filter_name, obs_date):
#     """
#     Save individual effective area plots.
#     """
#     filter_dir = OUTPUT_DIR / filter_name
#     filter_dir.mkdir(parents=True, exist_ok=True)
#     output_file = filter_dir / f"{filter_name}_{obs_date.replace(':', '').replace(' ', '_')}.png"
    
#     mismatch_count = len(failed_indices)

#     plt.figure(figsize=(10, 6))
#     plt.plot(wavelength, IDL_area, label="IDL Effective Area", color="blue", lw=2)
#     plt.plot(wavelength, XRTpy_area, label="XRTpy Effective Area", color="green", linestyle="--", lw=2)
#     plt.scatter(
#         wavelength[failed_indices],
#         IDL_area[failed_indices],
#         color="red",
#         label=f"Failed Points: {mismatch_count}",
#         zorder=5,
#     )
#     plt.xlabel("Wavelength (Å)", fontsize=14)
#     plt.ylabel("Effective Area (cm²)", fontsize=14)
#     plt.title(f"Effective Area Comparison for {filter_name} on {obs_date}", fontsize=16)
#     plt.legend(fontsize=12)
#     plt.grid(True, linestyle="--", alpha=0.7)
#     plt.tight_layout()
#     plt.savefig(output_file, dpi=300)
#     plt.close()
#     print(f"Saved plot to {output_file}")


# def save_combined_filter_plots(filter_name, test_data):
#     """
#     Save a combined plot with subplots for all dates associated with a single filter.
#     """
#     # Create a subdirectory for the filter
#     filter_dir = OUTPUT_DIR / filter_name
#     filter_dir.mkdir(parents=True, exist_ok=True)

#     # Define the output filename
#     output_file = filter_dir / f"{filter_name}_combined.png"

#     # Determine subplot grid size
#     num_tests = len(test_data)
#     num_cols = 3
#     num_rows = (num_tests + num_cols - 1) // num_cols  # Ceiling division

#     # Create the figure and axes
#     fig, axes = plt.subplots(num_rows, num_cols, figsize=(15, 5 * num_rows), squeeze=False)

#     for i, (wavelength, IDL_area, XRTpy_area, failed_indices, obs_date) in enumerate(test_data):
#         row, col = divmod(i, num_cols)
#         ax = axes[row, col]

#         # Plot data
#         ax.plot(wavelength, IDL_area, label="IDL Effective Area", color="blue", lw=2)
#         ax.plot(wavelength, XRTpy_area, label="XRTpy Effective Area", color="green", linestyle="--", lw=2)
#         ax.scatter(
#             wavelength[failed_indices],
#             IDL_area[failed_indices],
#             color="red",
#             label=f"Failed Points: {len(failed_indices)}",
#             zorder=5,
#         )

#         # Customize subplot
#         ax.set_title(f"{filter_name} - {obs_date}", fontsize=12)
#         ax.set_xlabel("Wavelength (Å)", fontsize=10)
#         ax.set_ylabel("Effective Area (cm²)", fontsize=10)
#         ax.grid(True, linestyle="--", alpha=0.7)
#         ax.legend(fontsize=8)

#     # Remove empty subplots
#     for j in range(i + 1, num_rows * num_cols):
#         fig.delaxes(axes[j // num_cols, j % num_cols])

#     # Adjust layout and save
#     plt.tight_layout()
#     plt.savefig(output_file, dpi=300)
#     plt.close()
#     print(f"Saved combined plot to {output_file}")


# @pytest.fixture(scope="session", autouse=True)
# def create_combined_plots():
#     """
#     Fixture to create combined plots for each filter after all tests are run.
#     """
#     yield  # Wait for all tests to complete

#     for filter_name, test_data in filter_test_data.items():
#         save_combined_filter_plots(filter_name, test_data)
########################################################################################################################
########################################################################################################################
########################################################################################################################


# import matplotlib.pyplot as plt

# @pytest.mark.parametrize("filename", get_IDL_data_files())
# def test_effective_area_compare_idl(filename):
#     print(f"\n\nTesting file: {filename}\n")
    
#     # Read the filter name and observation date from the file
#     with filename.open() as f:
#         filter_name = f.readline().split()[1]
#         filter_obs_date = " ".join(f.readline().split()[1:])
#         print(f"Filter name: {filter_name}, Observation date: {filter_obs_date}")  # Debugging output
    
#     # Correct non-standard date format
#     filter_obs_date = filter_obs_date.replace("Sept", "Sep")
    
#     # Load IDL data from the file
#     IDL_data = np.loadtxt(filename, skiprows=3)
#     IDL_wavelength = IDL_data[:, 0] * u.AA
#     IDL_effective_area = IDL_data[:, 1] * u.cm**2
    
#     # Compute effective area using XRTpy
#     instance = EffectiveAreaFundamental(filter_name, filter_obs_date)
#     actual_effective_area = instance.effective_area()
    
#     # Interpolate XRTpy effective area onto the IDL wavelength grid
#     XRTpy_effective_area_interp = np.interp(
#         IDL_wavelength.value,  # Target grid (IDL wavelengths)
#         instance.wavelength.value,  # Source grid (XRTpy wavelengths)
#         actual_effective_area.value  # Data to interpolate (XRTpy effective area)
#     )
    
#     # Make both arrays dimensionless for comparison
#     IDL_effective_area_unitless = IDL_effective_area.value
#     XRTpy_effective_area_unitless = XRTpy_effective_area_interp
    
#     # Compare effective areas using relative tolerance (unitless comparison)
#     rtol = 1e-4
#     differences = np.abs(XRTpy_effective_area_unitless - IDL_effective_area_unitless)
#     max_diff = np.max(differences)
#     failed_indices = np.where(differences > rtol * np.abs(IDL_effective_area_unitless))[0]
    
#     if failed_indices.size > 0:
#         print(f"Test failed for filter {filter_name} on {filter_obs_date}")
#         print(f"Max difference: {max_diff}")
        
#         # Plot the results
#         plot_effective_area_comparison(
#             IDL_wavelength.value,
#             IDL_effective_area_unitless,
#             XRTpy_effective_area_unitless,
#             failed_indices,
#             filter_name,
#             filter_obs_date,
#         )
        
#         # Fail the test
#         assert False, f"Effective areas differ for filter {filter_name} on {filter_obs_date}"


# def plot_effective_area_comparison(wavelength, IDL_area, XRTpy_area, failed_indices, filter_name, obs_date):
#     """
#     Plots the effective area comparison between XRTpy and IDL.
#     Marks failed points in red.
#     """
#     plt.figure(figsize=(10, 6))
#     plt.plot(wavelength, IDL_area, label="IDL Effective Area", color="blue", lw=2)
#     plt.plot(wavelength, XRTpy_area, label="XRTpy Effective Area", color="green", linestyle="--", lw=2)
#     plt.scatter(
#         wavelength[failed_indices],
#         IDL_area[failed_indices],
#         color="red",
#         label="Failed Points",
#         zorder=5,
#     )
#     plt.xlabel("Wavelength (Å)", fontsize=14)
#     plt.ylabel("Effective Area (cm²)", fontsize=14)
#     plt.title(f"Effective Area Comparison for {filter_name} on {obs_date}", fontsize=16)
#     plt.legend(fontsize=12)
#     plt.grid(True, linestyle="--", alpha=0.7)
#     plt.tight_layout()
#     plt.show()



########################################################################################################################
########################################################################################################################
########################################################################################################################
# Original 
## NOTE: This is marked as xfail because the IDL results that this test compares against
## are incorrect due to the use of quadratic interpolation in the contamination curves
## which leads to ringing near the edges in the contamination curve.
## See https://github.com/HinodeXRT/xrtpy/pull/284#issuecomment-2334503108
# @pytest.mark.xfail
# @pytest.mark.parametrize("filename", get_IDL_data_files())
# def test_effective_area_compare_idl(filename):
#     with Path.open(filename) as f:
#         filter_name = f.readline().split()[1]
#         filter_obs_date = " ".join(f.readline().split()[1:])
#     # NOTE: Annoyingly the date strings use "Sept" instead of "Sep" for "September"
#     filter_obs_date = filter_obs_date.replace("Sept", "Sep")
#     IDL_data = np.loadtxt(filename, skiprows=3)
#     IDL_wavelength = IDL_data[:, 0] * u.AA
#     IDL_effective_area = IDL_data[:, 1] * u.cm**2
#     instance = EffectiveAreaFundamental(filter_name, filter_obs_date)
#     actual_effective_area = instance.effective_area()
#     IDL_effective_area = np.interp(
#         instance.wavelength, IDL_wavelength, IDL_effective_area
#     )
#     assert u.allclose(
#         actual_effective_area,
#         IDL_effective_area,
#         rtol=1e-4, #Moderate Precision
#     )