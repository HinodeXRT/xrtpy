from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from utils_case_io import (
    read_mc_intensities_csv,
    run_dem_for_mc_csv,
    load_idl_dem_sav,
    case_dir,
)

#NOTE-User will need Python 3.11 to run

# CASE_DIR = case_dir("case_20210720_1604")
# csv_path = CASE_DIR / "mc_intensities_20210720_1604_IDL.csv"
# idl = load_idl_dem_sav(CASE_DIR / "xrt_dem_output_20210720_1604_MCITER100.sav")
# observation_date = "2021-07-20T16:04"

CASE_DIR = case_dir("case_20260107_124503")
csv_path = CASE_DIR / "mc_intensities_20260107_124503_IDL.csv"
idl = load_idl_dem_sav(CASE_DIR / "xrt_dem_output_20260107_124503_MCITER100.sav")
observation_date = "2026-01-07T12:45:03"


mc = read_mc_intensities_csv(csv_path)


print(mc.filters)
print(mc.mc_intensities.shape)  # (N, n_filters)
print(mc.df.head())



out = run_dem_for_mc_csv(
    csv_path=csv_path,
    observation_date=observation_date,
    intensity_errors=None,  # keep None to match IDL default behavior
)

#import pdb; pdb.set_trace()

print("filters:", out.filters)
print("logT shape:", out.logT.shape)
print("dem_runs shape:", out.dem_runs.shape)
print("modeled_runs shape:", out.modeled_runs.shape)
print("chisq_runs shape:", out.chisq_runs.shape)

#import pdb; pdb.set_trace()

print("\n\n-New Section-\n")
print("IDL logT:", idl.logT.shape)
print("IDL dem_base:", idl.dem_base.shape)
print("IDL dem_runs:", idl.dem_runs.shape, "n_runs=", idl.n_runs)


# import numpy as np
# import matplotlib.pyplot as plt
# from pathlib import Path

# def _log10_dem(dem: np.ndarray, floor: float = 1e-99) -> np.ndarray:
#     return np.log10(np.maximum(dem, floor))

# # -------------------------
# # Choose which 10 runs to plot
# # - Here: first 10 MC runs (Run 1..10)
# #   Your arrays are:
# #     idl.dem_runs: (100, 26)  -> MC runs only
# #     out.dem_runs: (100, 26)  -> MC runs only
# # So run number i corresponds to index i-1 in these arrays.
# # -------------------------
# runs_to_plot = list(range(0, 99))  # 1..10

# # Sanity checks
# assert idl.logT.shape == out.logT.shape
# assert np.allclose(idl.logT, out.logT, atol=1e-8), "IDL and Python logT grids differ!"
# logT = out.logT

# # Extract the 10 runs from each side
# idl_10 = np.stack([idl.dem_runs[i - 1] for i in runs_to_plot], axis=0)  # (10, 26)
# py_10  = np.stack([out.dem_runs[i - 1] for i in runs_to_plot], axis=0)  # (10, 26)

# log_idl_10 = _log10_dem(idl_10)
# log_py_10  = _log10_dem(py_10)

# # Global y-limits across all 10 overlays (same scale for every plot)
# ymin = np.min([log_idl_10.min(), log_py_10.min()])
# ymax = np.max([log_idl_10.max(), log_py_10.max()])

# # Optional padding for readability
# pad = 0.05 * (ymax - ymin) if ymax > ymin else 0.5
# ymin -= pad
# ymax += pad

# # Output folder
# plots_dir = CASE_DIR / "plots_overlay_first10"
# plots_dir.mkdir(parents=True, exist_ok=True)

# for k, run in enumerate(runs_to_plot):
#     fig, ax = plt.subplots(figsize=(8, 5))

#     ax.plot(logT, log_idl_10[k], label=f"IDL run {run}", linewidth=2)
#     ax.plot(logT, log_py_10[k],  label=f"XRTpy run {run}", linewidth=2, linestyle="--")

#     ax.set_title(f"DEM Overlay (Run {run})")
#     ax.set_xlabel("log T (K)")
#     ax.set_ylabel("log10(DEM)")
#     ax.set_ylim(ymin, ymax)
#     ax.grid(True, alpha=0.3)
#     ax.legend()

#     out_png = plots_dir / f"dem_overlay_run_{run:03d}.png"
#     fig.tight_layout()
#     fig.savefig(out_png, dpi=200)
#     plt.close(fig)


# print(f"\n\nSaved overlay plots to: {plots_dir}\n\n")


# # -----------------------------------------------------------------------------
# # SECTION 2: Make an MP4 movie (ffmpeg)
# # -----------------------------------------------------------------------------
# # If Matplotlib can't find ffmpeg automatically, uncomment and set path:
# # import matplotlib as mpl
# # mpl.rcParams["animation.ffmpeg_path"] = "/opt/homebrew/bin/ffmpeg"

# fig, ax = plt.subplots(figsize=(8, 5))

# line_idl, = ax.plot([], [], linewidth=2, label="IDL")
# line_py,  = ax.plot([], [], linewidth=2, linestyle="--", label="XRTpy")

# ax.set_xlim(float(logT.min()), float(logT.max()))
# ax.set_ylim(ymin, ymax)
# ax.set_xlabel("log T (K)")
# ax.set_ylabel("log10(DEM)")
# ax.grid(True, alpha=0.3)
# ax.legend()

# def update(frame: int):
#     run_idx = runs_to_plot[frame]
#     ax.set_title(f"DEM Overlay (Run {run_idx + 1})")

#     line_idl.set_data(logT, log_idl[frame])
#     line_py.set_data(logT, log_py[frame])

#     return line_idl, line_py

# ani = animation.FuncAnimation(
#     fig,
#     update,
#     frames=len(runs_to_plot),
#     blit=True,
# )

# movie_path = CASE_DIR / "dem_overlay.mp4"

# writer = animation.FFMpegWriter(
#     fps=2,  # change speed here
#     metadata={"artist": "xrtpy"},
#     bitrate=1800,
# )

# ani.save(movie_path, writer=writer, dpi=200)
# plt.close(fig)

# print(f"Movie saved to: {movie_path}")



def _log10_dem(dem: np.ndarray, floor: float = 1e-99) -> np.ndarray:
    return np.log10(np.maximum(dem, floor))

# -------------------------
# Choose runs to plot
# Use 1-based run numbers for filenames/titles, but index arrays with run-1
# -------------------------
runs_to_plot = list(range(1, 101))  # 1..10 (matches your "first 10" folder)
# If you want all 100:
# runs_to_plot = list(range(1, 101))  # 1..100

# Sanity checks
assert idl.logT.shape == out.logT.shape
assert np.allclose(idl.logT, out.logT, atol=1e-8), "IDL and Python logT grids differ!"
logT = out.logT

# Extract selected runs from each side
idl_sel = np.stack([idl.dem_runs[r - 1] for r in runs_to_plot], axis=0)  # (n, 26)
py_sel  = np.stack([out.dem_runs[r - 1] for r in runs_to_plot], axis=0)  # (n, 26)

log_idl_10 = _log10_dem(idl_sel)
log_py_10  = _log10_dem(py_sel)

# Global y-limits across all selected overlays (same scale)
ymin = float(min(log_idl_10.min(), log_py_10.min()))
ymax = float(max(log_idl_10.max(), log_py_10.max()))
pad = 0.05 * (ymax - ymin) if ymax > ymin else 0.5
ymin -= pad
ymax += pad

# -------------------------
# SECTION 1: Save overlay PNGs
# -------------------------
plots_dir = CASE_DIR / "plots_overlay"
plots_dir.mkdir(parents=True, exist_ok=True)

for k, run in enumerate(runs_to_plot):
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(logT, log_idl_10[k], label=f"IDL run {run}", linewidth=2)
    ax.plot(logT, log_py_10[k],  label=f"XRTpy run {run}", linewidth=2, linestyle="--")

    ax.set_title(f"DEM Overlay (Run {run})")
    ax.set_xlabel("log T (K)")
    ax.set_ylabel("log10(DEM)")
    ax.set_ylim(ymin, ymax)
    ax.grid(True, alpha=0.3)
    ax.legend()

    out_png = plots_dir / f"dem_overlay_run_{run:03d}.png"
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

print(f"Saved {len(runs_to_plot)} overlay plots to: {plots_dir}")

# -------------------------
# SECTION 2: Make MP4 movie via ffmpeg
# -------------------------
# If Matplotlib can't find ffmpeg automatically, uncomment:
# import matplotlib as mpl
# mpl.rcParams["animation.ffmpeg_path"] = "/opt/homebrew/bin/ffmpeg"

fig, ax = plt.subplots(figsize=(8, 5))

(line_idl,) = ax.plot([], [], linewidth=2, label="IDL")
(line_py,)  = ax.plot([], [], linewidth=2, linestyle="--", label="XRTpy")

ax.set_xlim(float(logT.min()), float(logT.max()))
ax.set_ylim(ymin, ymax)
ax.set_xlabel("log T (K)")
ax.set_ylabel("log10(DEM)")
ax.grid(True, alpha=0.3)
ax.legend()

def update(frame: int):
    run = runs_to_plot[frame]
    ax.set_title(f"DEM Overlay (Run {run})")
    line_idl.set_data(logT, log_idl_10[frame])
    line_py.set_data(logT,  log_py_10[frame])
    return (line_idl, line_py)

ani = animation.FuncAnimation(
    fig,
    update,
    frames=len(runs_to_plot),
    blit=False,  # safer (avoids the blit init crash you saw)
)

movie_path = plots_dir / "dem_overlay.mp4"

writer = animation.FFMpegWriter(
    fps=6.5,
    metadata={"artist": "xrtpy"},
    bitrate=1800,
)

ani.save(movie_path, writer=writer, dpi=200)
plt.close(fig)

print(f"Movie saved to: {movie_path}")

######

# -----------------------------------------------------------------------------
# Professional PNGs (2 panels) that MATCH the movie style:
#  - Top: log10(DEM) overlay (IDL vs XRTpy)
#  - Bottom: Δ(log10 DEM) at major logT points with ±1σ error bars (across runs)
#     Δ(log10 DEM) = log10(DEM_XRTpy) - log10(DEM_IDL)
# -----------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

def _log10_dem(dem: np.ndarray, floor: float = 1e-99) -> np.ndarray:
    return np.log10(np.maximum(dem, floor))

# runs_to_plot should be 1-based run numbers, e.g.:
# runs_to_plot = list(range(1, 11))   # first 10
# runs_to_plot = list(range(1, 101))  # all 100

# logT should already exist, and log_idl_10/log_py_10 should correspond to runs_to_plot order
# If you need to rebuild them:
# idl_sel = np.stack([idl.dem_runs[r - 1] for r in runs_to_plot], axis=0)
# py_sel  = np.stack([out.dem_runs[r - 1] for r in runs_to_plot], axis=0)
# log_idl_10 = _log10_dem(idl_sel)
# log_py_10  = _log10_dem(py_sel)

# -------------------------
# Major logT points (5.0, 5.5, ... 8.0) and nearest indices
# -------------------------
major_logT = np.arange(5.0, 8.0 + 1e-9, 0.5)
major_idx = np.array([int(np.argmin(np.abs(logT - t))) for t in major_logT])
major_x = logT[major_idx]  # actual grid values used

# -------------------------
# Δ(log10 DEM) across ALL selected runs
# -------------------------
delta_all = log_py_10 - log_idl_10                     # (n_runs_selected, 26)
sigma_delta = np.std(delta_all, axis=0)                # (26,)
sigma_major = sigma_delta[major_idx]                   # (n_major,)

# y-limits for bottom panel based on all runs at major points (stable scaling)
dmin = float(np.min(delta_all[:, major_idx]))
dmax = float(np.max(delta_all[:, major_idx]))
dpad = 0.15 * (dmax - dmin) if dmax > dmin else 0.5
dmin -= dpad
dmax += dpad

# Optional "agreement band" shading (dex). Set to None to disable.
delta_band = None          # e.g. 0.10 for ±0.10 dex band, or None to disable

# -------------------------
# SECTION 1: Save professional PNGs (2 panels)
# -------------------------
plots_dir = CASE_DIR / "plots_overlay_w_diff"
plots_dir.mkdir(parents=True, exist_ok=True)

for k, run in enumerate(runs_to_plot):
    # run is 1-based; k is 0-based index into log_idl_10/log_py_10/delta_all
    chisq = float(out.chisq_runs[run - 1]) if hasattr(out, "chisq_runs") else np.nan

    fig, (ax1, ax2) = plt.subplots(
        2, 1,
        figsize=(8.5, 7.0),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
        constrained_layout=True,
    )

    # ---- Top: DEM overlay ----
    ax1.plot(logT, log_idl_10[k], linewidth=2.2, label="IDL")
    ax1.plot(logT, log_py_10[k],  linewidth=2.2, linestyle="--", label="XRTpy")
    ax1.set_ylabel("log10(DEM)")
    ax1.set_ylim(ymin, ymax)
    ax1.grid(True, alpha=0.25)
    ax1.legend(loc="best", frameon=True)

    # Annotate run + chisq
    ax1.text(
        0.02, 0.95,
        f"Run: {run}\n$\\chi^2$: {chisq:.3g}",
        transform=ax1.transAxes,
        va="top", ha="left",
        bbox=dict(boxstyle="round", alpha=0.15),
    )

    # ---- Bottom: Δ(log10 DEM) at major points with ±1σ error bars ----
    ax2.axhline(0.0, linewidth=1.3, alpha=0.8)

    y_major = delta_all[k, major_idx]  # Δ at major points for this run

    ax2.errorbar(
        major_x,
        y_major,
        yerr=sigma_major,
        fmt="o",
        capsize=3,
        elinewidth=1.0,
        markersize=4.5,
        color="black",
    )

    if delta_band is not None:
        ax2.axhspan(-delta_band, +delta_band, alpha=0.12)

    ax2.set_xlabel("log T (K)")
    ax2.set_ylabel("Δ log10(DEM)\n(XRTpy − IDL)")
    ax2.set_ylim(dmin, dmax)
    ax2.grid(True, alpha=0.25)

    fig.suptitle("DEM Comparison: IDL vs XRTpy", y=1.02)

    out_png = plots_dir / f"dem_compare_run_{run:03d}.png"
    fig.savefig(out_png, dpi=250)
    plt.close(fig)

print(f"Saved {len(runs_to_plot)} professional plots to: {plots_dir}")
# # -------------------------
# # SECTION 2: Movie (MP4) with same professional layout
# # -------------------------
# fig, (ax1, ax2) = plt.subplots(
#     2, 1, figsize=(8.5, 7.0),
#     sharex=True,
#     gridspec_kw={"height_ratios": [3, 1]},
#     constrained_layout=True,
# )

# # Pre-create artists (faster + cleaner)
# (line_idl,) = ax1.plot([], [], linewidth=2.2, label="IDL")
# (line_py,)  = ax1.plot([], [], linewidth=2.2, linestyle="--", label="XRTpy")
# ax1.set_ylabel("log10(DEM)")
# ax1.set_ylim(ymin, ymax)
# ax1.grid(True, alpha=0.25)
# ax1.legend(loc="best", frameon=True)

# # annotation text in top panel
# anno = ax1.text(
#     0.02, 0.95, "",
#     transform=ax1.transAxes,
#     va="top", ha="left",
#     bbox=dict(boxstyle="round", alpha=0.15),
# )

# ax2.axhline(0.0, linewidth=1.5, alpha=0.7)
# (line_delta,) = ax2.plot([], [], linewidth=2.0)
# ax2.set_xlabel("log10(T [K])")
# ax2.set_ylabel("Δ dex")
# ax2.set_ylim(dmin, dmax)
# ax2.grid(True, alpha=0.25)

# if delta_band is not None:
#     ax2.axhspan(-delta_band, +delta_band, alpha=0.12)

# ax2.set_xlim(float(logT.min()), float(logT.max()))
# fig.suptitle("DEM Comparison: IDL vs XRTpy", y=1.02)

# def update(frame: int):
#     run = runs_to_plot[frame]  # 1-based
#     chisq = float(out.chisq_runs[run - 1]) if hasattr(out, "chisq_runs") else np.nan

#     line_idl.set_data(logT, log_idl_10[frame])
#     line_py.set_data(logT,  log_py_10[frame])
#     line_delta.set_data(logT, delta[frame])

#     anno.set_text(f"Run: {run}\n$\\chi^2$: {chisq:.3g}")
#     return (line_idl, line_py, line_delta, anno)

# ani = animation.FuncAnimation(fig, update, frames=len(runs_to_plot), blit=False)

# # Speed control
# fps = 15  # increase for faster video
# movie_path = plots_dir / "dem_compare_professional_v2.mp4"
# writer = animation.FFMpegWriter(fps=fps, metadata={"artist": "xrtpy"}, bitrate=2200)

# ani.save(movie_path, writer=writer, dpi=250)
# plt.close(fig)

# print(f"Movie saved to: {movie_path}")
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import LineCollection

def _log10_dem(dem: np.ndarray, floor: float = 1e-99) -> np.ndarray:
    return np.log10(np.maximum(dem, floor))

# -------------------------
# Inputs assumed to exist:
#   logT (nT,)
#   log_idl_10 (n_runs, nT)
#   log_py_10  (n_runs, nT)
#   runs_to_plot (list of run numbers, 1-based)
#   ymin, ymax
#   plots_dir (Path)
# -------------------------

# Major points: only keep ones inside the grid
major_logT = np.arange(5.0, 8.0 + 1e-9, 0.5)
major_logT = major_logT[(major_logT >= logT.min() - 1e-9) & (major_logT <= logT.max() + 1e-9)]

major_idx = np.array([int(np.argmin(np.abs(logT - t))) for t in major_logT])
major_x = logT[major_idx]

# Delta across frames
delta_all = log_py_10 - log_idl_10                  # (n_frames, nT)
sigma_delta = np.std(delta_all, axis=0)             # (nT,)
sigma_major = sigma_delta[major_idx]                # (n_major,)

# Stable y-range for bottom panel
dmin = float(np.min(delta_all[:, major_idx]))
dmax = float(np.max(delta_all[:, major_idx]))
pad = 0.15 * (dmax - dmin) if dmax > dmin else 0.5
ax_bot_ylim = (dmin - pad, dmax + pad)

# -------------------------
# Figure + axes
# -------------------------
fig = plt.figure(figsize=(9, 7))
gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[3, 1], hspace=0.12)

ax_top = fig.add_subplot(gs[0, 0])
ax_bot = fig.add_subplot(gs[1, 0], sharex=ax_top)

# Top lines
(line_idl,) = ax_top.plot([], [], linewidth=2, label="IDL")
(line_py,)  = ax_top.plot([], [], linewidth=2, linestyle="--", label="XRTpy")

ax_top.set_xlim(float(logT.min()), float(logT.max()))
ax_top.set_ylim(ymin, ymax)
ax_top.set_ylabel("log10(DEM)")
ax_top.grid(True, alpha=0.3)
ax_top.legend(loc="best")

# Bottom baseline
ax_bot.axhline(0.0, linewidth=1, alpha=0.8)
ax_bot.set_xlabel("log T (K)")
ax_bot.set_ylabel("Δ log10(DEM)\n(XRTpy - IDL)")
ax_bot.grid(True, alpha=0.3)
ax_bot.set_ylim(*ax_bot_ylim)

plt.setp(ax_top.get_xticklabels(), visible=False)

# -------------------------
# Build animated "errorbar" artists manually
# -------------------------
# Markers
(points_line,) = ax_bot.plot(major_x, np.zeros_like(major_x), "o", markersize=5, color="black")

# Vertical error segments (LineCollection)
def make_vsegments(y_major):
    # segments shape: (n_major, 2, 2) => [ [ [x,ylo],[x,yhi] ], ...]
    ylo = y_major - sigma_major
    yhi = y_major + sigma_major
    segs = np.zeros((len(major_x), 2, 2), dtype=float)
    segs[:, 0, 0] = major_x
    segs[:, 1, 0] = major_x
    segs[:, 0, 1] = ylo
    segs[:, 1, 1] = yhi
    return segs

# Caps (small horizontal lines at top/bottom)
cap_halfwidth = 0.03  # in logT units (adjust if you want)
def make_capsegs(y_major):
    ylo = y_major - sigma_major
    yhi = y_major + sigma_major

    segs = []
    for x, yl, yh in zip(major_x, ylo, yhi):
        segs.append([[x - cap_halfwidth, yl], [x + cap_halfwidth, yl]])
        segs.append([[x - cap_halfwidth, yh], [x + cap_halfwidth, yh]])
    return np.array(segs, dtype=float)

vlines = LineCollection(make_vsegments(np.zeros_like(major_x)), colors="black", linewidths=1.0)
caps   = LineCollection(make_capsegs(np.zeros_like(major_x)), colors="black", linewidths=1.0)
ax_bot.add_collection(vlines)
ax_bot.add_collection(caps)

# -------------------------
# Update function
# -------------------------
def update(frame: int):
    run = runs_to_plot[frame]  # 1-based run label

    # Top
    ax_top.set_title(f"DEM Overlay (Run {run})")
    line_idl.set_data(logT, log_idl_10[frame])
    line_py.set_data(logT,  log_py_10[frame])

    # Bottom
    y_major = delta_all[frame, major_idx]
    points_line.set_data(major_x, y_major)

    vlines.set_segments(make_vsegments(y_major))
    caps.set_segments(make_capsegs(y_major))

    return (line_idl, line_py, points_line, vlines, caps)

ani = animation.FuncAnimation(
    fig,
    update,
    frames=len(runs_to_plot),
    blit=False,  # keep robust
)

movie_path = plots_dir / "dem_overlay_with_diff.mp4"
writer = animation.FFMpegWriter(
    fps=7.5,          # increase -> faster
    metadata={"artist": "xrtpy"},
    bitrate=1800,
)

ani.save(movie_path, writer=writer, dpi=200)
plt.close(fig)

print(f"Movie saved to: {movie_path}")
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation

# # -------------------------
# # Bottom-panel "major" logT points
# # -------------------------
# major_logT = np.arange(5.0, 8.0 + 1e-9, 0.5)  # 5.0, 5.5, ... 8.0

# # For each major temp, pick the nearest index in the logT grid
# major_idx = np.array([int(np.argmin(np.abs(logT - t))) for t in major_logT])
# major_x = logT[major_idx]  # actual grid values used (closest to requested)

# # Compute delta(logDEM) for all frames/runs we are animating
# # shape: (n_frames, 26)
# delta_all = log_py_10 - log_idl_10

# # Error bars: 1-sigma across frames at each temperature bin (professional context)
# # shape: (26,)
# sigma_delta = np.std(delta_all, axis=0)

# # Only at the major points:
# sigma_major = sigma_delta[major_idx]

# # -------------------------
# # Figure with two subplots
# # -------------------------
# fig = plt.figure(figsize=(9, 7))
# gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[3, 1], hspace=0.12)

# ax_top = fig.add_subplot(gs[0, 0])
# ax_bot = fig.add_subplot(gs[1, 0], sharex=ax_top)

# # Top lines
# (line_idl,) = ax_top.plot([], [], linewidth=2, label="IDL")
# (line_py,)  = ax_top.plot([], [], linewidth=2, linestyle="--", label="XRTpy")

# ax_top.set_xlim(float(logT.min()), float(logT.max()))
# ax_top.set_ylim(ymin, ymax)
# ax_top.set_ylabel("log10(DEM)")
# ax_top.grid(True, alpha=0.3)
# ax_top.legend(loc="best")

# # Bottom: delta points with error bars (we'll update their y-values each frame)
# # We create an "errorbar container" once and then update the data each frame.
# eb = ax_bot.errorbar(
#     major_x,
#     np.zeros_like(major_x),
#     yerr=sigma_major,
#     fmt="o",
#     capsize=3,
#     elinewidth=1,
#     color='Black'
# )

# ax_bot.axhline(0.0, linewidth=1, alpha=0.8)
# ax_bot.set_xlabel("log T (K)")
# ax_bot.set_ylabel("Δ log10(DEM)\n(XRTpy - IDL)")
# ax_bot.grid(True, alpha=0.3)

# # Set a reasonable fixed y-range for delta based on overall spread
# # (prevents the bottom axis from jumping around frame-to-frame)
# dmin = float(np.min(delta_all[:, major_idx]))
# dmax = float(np.max(delta_all[:, major_idx]))
# pad = 0.15 * (dmax - dmin) if dmax > dmin else 0.5
# ax_bot.set_ylim(dmin - pad, dmax + pad)

# # Prevent top subplot from repeating x tick labels
# plt.setp(ax_top.get_xticklabels(), visible=False)

# def update(frame: int):
#     run = runs_to_plot[frame]  # run number (1-based in your setup)

#     # Top panel
#     ax_top.set_title(f"DEM Overlay (Run {run})")
#     line_idl.set_data(logT, log_idl_10[frame])
#     line_py.set_data(logT,  log_py_10[frame])

#     # Bottom panel: delta at major temp points
#     y_major = delta_all[frame, major_idx]

#     # Update errorbar artist:
#     # eb.lines[0] is the marker/line for the data points
#     eb.lines[0].set_data(major_x, y_major)

#     # Return artists for blitting (we'll keep blit=False for robustness)
#     return (line_idl, line_py, eb.lines[0])

# ani = animation.FuncAnimation(
#     fig,
#     update,
#     frames=len(runs_to_plot),
#     blit=False,
# )

# movie_path = plots_dir / "dem_overlay_with_diff.mp4"

# writer = animation.FFMpegWriter(
#     fps=7.5,  # <-- make faster/slower here
#     metadata={"artist": "xrtpy"},
#     bitrate=1800,
# )

# ani.save(movie_path, writer=writer, dpi=200)
# plt.close(fig)

# print(f"Science-style movie saved to: {movie_path}")
######

# Single plot: overlay ALL XRTpy DEM runs on one figure (log10 DEM vs logT)
# Assumes you already have:
#   - out.dem_runs shape (n_runs, 26)  (here n_runs=100)
#   - out.logT shape (26,)

import numpy as np
import matplotlib.pyplot as plt

def _log10_dem(dem: np.ndarray, floor: float = 1e-99) -> np.ndarray:
    return np.log10(np.maximum(dem, floor))

logT = out.logT
dem_runs = out.dem_runs  # (100, 26)

log_dem_runs = _log10_dem(dem_runs)

fig, ax = plt.subplots(figsize=(9, 6))

# Plot every run (thin)
for i in range(log_dem_runs.shape[0]):
    ax.plot(logT, log_dem_runs[i], linewidth=1.0, alpha=0.25)


ax.set_title("XRTpy DEM Monte Carlo Runs (All Overlaid)")
ax.set_xlabel("log T (K)")
ax.set_ylabel("log10(DEM)")
ax.grid(True, alpha=0.3)
ax.legend(loc="best")

# Save
out_png = CASE_DIR / "xrtpy_dem_all_runs_overlay.png"
fig.tight_layout()
fig.savefig(out_png, dpi=250)
plt.close(fig)

print(f"Saved: {out_png}")
######
# Single plot: overlay ALL IDL DEM runs (log10 DEM vs logT)
# Assumes:
#   - idl.dem_runs shape (100, 26)
#   - idl.logT shape (26,)

import numpy as np
import matplotlib.pyplot as plt

def _log10_dem(dem: np.ndarray, floor: float = 1e-99) -> np.ndarray:
    return np.log10(np.maximum(dem, floor))

logT = idl.logT
dem_runs = idl.dem_runs  # (100, 26)

log_dem_runs = _log10_dem(dem_runs)

fig, ax = plt.subplots(figsize=(9, 6))

# Plot every IDL run
for i in range(log_dem_runs.shape[0]):
    ax.plot(logT, log_dem_runs[i], linewidth=1.0, alpha=0.35)

ax.set_title("IDL DEM Monte Carlo Runs (All Overlaid)")
ax.set_xlabel("log T (K)")
ax.set_ylabel("log10(DEM)")
ax.grid(True, alpha=0.3)

# Save
out_png = CASE_DIR / "idl_dem_all_runs_overlay.png"
fig.tight_layout()
fig.savefig(out_png, dpi=250)
plt.close(fig)

print(f"Saved: {out_png}")
######
import numpy as np
import matplotlib.pyplot as plt

def _log10_dem(dem: np.ndarray, floor: float = 1e-99) -> np.ndarray:
    return np.log10(np.maximum(dem, floor))

# Ensure grids match
assert np.allclose(idl.logT, out.logT, atol=1e-8), "logT grids differ!"
logT = idl.logT

log_idl = _log10_dem(idl.dem_runs)      # (100, 26)
log_xrt = _log10_dem(out.dem_runs)     # (100, 26)

fig, ax = plt.subplots(figsize=(9, 6))

# --- Plot IDL (orange) ---
for i in range(log_idl.shape[0]):
    ax.plot(logT, log_idl[i], color="orange", linewidth=1.0, alpha=0.35)

# --- Plot XRTpy (blue) ---
for i in range(log_xrt.shape[0]):
    ax.plot(logT, log_xrt[i], color="blue", linewidth=1.0, alpha=0.35)

ax.set_title("DEM Monte Carlo Runs: IDL (Orange) vs XRTpy (Blue)")
ax.set_xlabel("log T (K)")
ax.set_ylabel("log10(DEM)")
ax.grid(True, alpha=0.3)

# Create manual legend entries (so we don’t get 200 legend lines)
from matplotlib.lines import Line2D
legend_lines = [
    Line2D([0], [0], color="orange", lw=2, label="IDL"),
    Line2D([0], [0], color="blue", lw=2, label="XRTpy"),
]
ax.legend(handles=legend_lines, loc="best")

fig.tight_layout()

out_png = CASE_DIR / "idl_vs_xrtpy_all_runs_overlay.png"
fig.savefig(out_png, dpi=250)
plt.close(fig)

print(f"Saved: {out_png}")
######

import numpy as np

print("\n--- INPUTS USED TO FEED XRTpy DEM ---")
print("observation_date:", observation_date)
print("csv_path:", csv_path)
print("filters (order):", mc.filters)
print("mc.mc_intensities shape:", mc.mc_intensities.shape)

# What is the "base" / "nominal" row?
# In many pipelines row 0 is the nominal intensities, and rows 1..N-1 are MC perturbed.
base_idx = 0

base_intensities = np.asarray(mc.mc_intensities[base_idx], dtype=float)
print("\nBase (row 0) intensities used for XRTpy:")
for f, val in zip(mc.filters, base_intensities):
    print(f"  {f:>12s}: {val:.6f} DN/s")

print("\nFirst 3 MC rows (sanity):")
for i in range(3):
    row = np.asarray(mc.mc_intensities[i], dtype=float)
    print(f"  row {i}: " + ", ".join([f"{v:.3f}" for v in row]))



import numpy as np
import matplotlib.pyplot as plt

def _log10_dem(dem: np.ndarray, floor: float = 1e-99) -> np.ndarray:
    return np.log10(np.maximum(dem, floor))

# Ensure temperature grids match
assert np.allclose(idl.logT, out.logT, atol=1e-8), "logT grids differ!"
logT = idl.logT

# --- IDL initial DEM ---
idl_initial = idl.dem_base  # from .sav file
log_idl_initial = _log10_dem(idl_initial)

# --- XRTpy initial DEM ---
# If your class defines dem_base, use it.
# Otherwise fallback to first run.
if hasattr(out, "dem_base"):
    xrt_initial = out.dem_base
else:
    xrt_initial = out.dem_runs[0]

log_xrt_initial = _log10_dem(xrt_initial)

# -------------------------
# Plot
# -------------------------
fig, ax = plt.subplots(figsize=(8.5, 6))

ax.plot(logT, log_idl_initial, color="orange", linewidth=2.5, label="IDL (initial)")
ax.plot(logT, log_xrt_initial, color="blue", linewidth=2.5, linestyle="--",
        label="XRTpy (initial)")

ax.set_title("Initial DEM Comparison: IDL vs XRTpy")
ax.set_xlabel("log T (K)")
ax.set_ylabel("log10(DEM)")
ax.grid(True, alpha=0.3)
ax.legend(loc="best")

fig.tight_layout()

out_png = CASE_DIR / "initial_dem_idl_vs_xrtpy.png"
fig.savefig(out_png, dpi=300)
plt.close(fig)

print(f"Saved: {out_png}")
