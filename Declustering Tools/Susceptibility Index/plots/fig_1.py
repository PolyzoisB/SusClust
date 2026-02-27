import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter


def histogram2d_density(x, y, bins, xlim, ylim):
    """
    Returns a 2D density estimate using a 2D histogram.

    x, y: 1D arrays
    bins: (nx, ny)
    xlim, ylim: (min, max)
    """
    nx, ny = bins
    H, xedges, yedges = np.histogram2d(
        x, y, bins=[nx, ny], range=[xlim, ylim]
    )
    # Convert counts to density (per unit area)
    dx = (xedges[-1] - xedges[0]) / nx
    dy = (yedges[-1] - yedges[0]) / ny
    density = H / (np.sum(H) * dx * dy) if np.sum(H) > 0 else H

    # Return centers (for pcolormesh we will use edges)
    return density, xedges, yedges


# -------------------- Load Susceptibility Index --------------------
df_si = pd.read_csv("../results/susc_index.txt", sep=r"\s+", header=None, engine="python")
df_si.columns = ["K", "SI", "NN", "Threshold", "ThresholdIndex", "Nlinks"]
df_si = df_si.astype(float)

K = df_si["K"].values
SI = df_si["SI"].values
NN = df_si["NN"].values
thr = df_si["Threshold"].values

# -------------------- Load final results --------------------
df_final = pd.read_csv("../results/final_results.txt", sep=r"\s+", header=None, engine="python")
# The Fortran writes 10 values; keep generic column names
df_final = df_final.astype(float)
K_best = int(df_final.iloc[0, 0])  # shifted-minimum clusters
thr_best = float(df_final.iloc[0, 2])

# -------------------- Load rescaled distances --------------------
df_dist = pd.read_csv("../results/rescaled_distances.txt", sep=r"\s+", header=None, engine="python")
df_dist.columns = ["T", "R"]
df_dist = df_dist.astype(float)

TN = df_dist["T"].values
RN = df_dist["R"].values

log_TN = np.log10(TN)
log_RN = np.log10(RN)

# 2D density on log-log plane
xlim = (-6, 7)
ylim = (-9, 1)
density, xedges, yedges = histogram2d_density(log_TN, log_RN, bins=(200, 200), xlim=xlim, ylim=ylim)

# Smooth density a bit (as in original)
density_s = gaussian_filter(density.T, sigma=1.5)  # transpose so axes match pcolormesh usage later

# Mask lowest 5% cumulative density (original behavior)
flat = density_s[np.isfinite(density_s)].ravel()
if flat.size > 0:
    s = np.sort(flat)
    c = np.cumsum(s)
    cutoff = s[np.searchsorted(c, 0.05 * c[-1])]
    density_s[density_s < cutoff] = np.nan

# Lines for threshold (kept consistent with your original approach)
Ta_centers = 0.5 * (xedges[:-1] + xedges[1:])
a_best = np.log10(1 / thr_best)
line_best = a_best * np.ones_like(Ta_centers) - Ta_centers

# Compute plot bounds from finite density
finite = np.isfinite(density_s)
if np.any(finite):
    x_min = 10 ** xedges[0]
    x_max = 10 ** xedges[-1]
    y_min = 10 ** yedges[0]
    y_max = 10 ** yedges[-1]
else:
    x_min, x_max = 1e-6, 1e7
    y_min, y_max = 1e-9, 1e1

# -------------------- Plot --------------------
plt.rcParams.update({"font.size": 14})
fig = plt.figure(figsize=(21, 8), dpi=400)
fig.patch.set_facecolor("white")

# Panel A: rescaled distances density
ax1 = plt.subplot(1, 2, 1)
# pcolormesh expects edges
X = 10 ** xedges
Y = 10 ** yedges
h = ax1.pcolormesh(X, Y, density_s, shading="auto")
ax1.plot(10 ** Ta_centers, 10 ** line_best, linewidth=3, linestyle="--", color="m",
         label=rf"$K_{{best}}={K_best}$")

plt.colorbar(h, ax=ax1, label="Density")
ax1.set_title("Rescaled distances")
ax1.set_xlabel("Rescaled time, T")
ax1.set_ylabel("Rescaled distance, R")
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlim([x_min, x_max])
ax1.set_ylim([y_min, y_max])
ax1.legend()

# Panel B: SI + NN histogram
ax2 = plt.subplot(1, 2, 2)
SI_norm = SI / np.nansum(SI) if np.nansum(SI) > 0 else SI
NN_pos = NN[NN > 0]
K_pos = K[NN > 0]
NN_norm = NN_pos / np.sum(NN_pos) if NN_pos.size > 0 else NN_pos

ax2.plot(K, SI_norm, linewidth=1.5, marker="o", color="b", label="SI / sum(SI)")
ax2.bar(K_pos, NN_norm, width=1, align="edge",
        facecolor="green", edgecolor="green", linewidth=1, alpha=0.5, label="H(K) / sum(H)")

ax2.axvline(x=K_best, linewidth=3, color="m", linestyle="--", label=rf"$K_{{best}}$")

ax2.set_title("NN Histogram")
ax2.set_xlabel("K")
ax2.legend()

plt.savefig("Fig_1_SI.png", bbox_inches="tight", dpi=400)
