import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter


def histogram2d_density(x, y, bins, xedges, yedges):
    H, _, _ = np.histogram2d(x, y, bins=[xedges, yedges])

    # Convert counts to density per unit area in (x,y) space
    dx = np.diff(xedges)
    dy = np.diff(yedges)
    area = dx[:, None] * dy[None, :]
    total = np.sum(H)
    density = (H / total) / area if total > 0 else H
    return density


# -------------------- Load Susceptibility Index --------------------
df_si = pd.read_csv("../results/susc_index.txt", sep=r"\s+", header=None, engine="python")
df_si.columns = ["K", "SI", "NN", "Threshold", "ThresholdIndex", "Nlinks"]
df_si = df_si.astype(float)

K = df_si["K"].values
SI = df_si["SI"].values
NN = df_si["NN"].values

# -------------------- Load final results --------------------
df_final = pd.read_csv("../results/final_results.txt", sep=r"\s+", header=None, engine="python").astype(float)
K_best = int(df_final.iloc[0, 0])
# Similarity threshold in sec*km (not rescaled)
thr_best = float(df_final.iloc[0, 2])

# -------------------- Load rescaled distances --------------------
df_dist = pd.read_csv("../results/rescaled_distances.txt", sep=r"\s+", header=None, engine="python")
df_dist.columns = ["T", "R"]
df_dist = df_dist.astype(float)

TN = df_dist["T"].values
RN = df_dist["R"].values

# log10-space
log_TN = np.log10(TN)
log_RN = np.log10(RN)

# ---- Data-driven bounds (avoids "out of bounds" plots) ----
# Use robust percentiles to ignore extreme outliers
xlo, xhi = np.nanpercentile(log_TN, [0.5, 99.5])
ylo, yhi = np.nanpercentile(log_RN, [0.5, 99.5])

# Add a small padding
padx = 0.05 * (xhi - xlo) if xhi > xlo else 1.0
pady = 0.05 * (yhi - ylo) if yhi > ylo else 1.0
xlo, xhi = xlo - padx, xhi + padx
ylo, yhi = ylo - pady, yhi + pady

nx = ny = 200
xedges = np.linspace(xlo, xhi, nx + 1)
yedges = np.linspace(ylo, yhi, ny + 1)

density = histogram2d_density(log_TN, log_RN, bins=(nx, ny), xedges=xedges, yedges=yedges)

# smooth (transpose for plotting with pcolormesh later)
density_s = gaussian_filter(density.T, sigma=1.5)

# Mask lowest 5% cumulative density
flat = density_s[np.isfinite(density_s)].ravel()
if flat.size > 0:
    s = np.sort(flat)
    c = np.cumsum(s)
    cutoff = s[np.searchsorted(c, 0.001 * c[-1])]
    density_s[density_s < cutoff] = np.nan

# Threshold line on log-log plane
xcenters = 0.5 * (xedges[:-1] + xedges[1:])
a_best = np.log10(1 / thr_best)
line_best = a_best * np.ones_like(xcenters) - xcenters
# Convert axes edges back to linear for plotting
X = 10 ** xedges
Y = 10 ** yedges

# -------------------- Plot --------------------
plt.rcParams.update({"font.size": 16})
fig = plt.figure(figsize=(21, 8), dpi=400)
fig.patch.set_facecolor("white")

# Panel A
ax1 = plt.subplot(1, 2, 1)
h = ax1.pcolormesh(X, Y, density_s, shading="auto")
ax1.plot(10 ** xcenters, 10 ** line_best, linewidth=3, linestyle="--", color="m",
         label=rf"$K_{{best}}={K_best}$")
plt.colorbar(h, ax=ax1, label="Density")

ax1.set_title("Rescaled distances", fontsize=14)
ax1.set_xlabel("Rescaled time, T (years)", fontsize=14)
ax1.set_ylabel("Rescaled distance, R (km)", fontsize=14)
ax1.set_xscale("log")
ax1.set_yscale("log")

# IMPORTANT: limits derived from the same edges used for pcolormesh
ax1.set_xlim([X[0], X[-1]])
ax1.set_ylim([Y[0], Y[-1]])
ax1.legend()
# --- Diagonal threshold line spanning the full axes ---
# Use current x-limits (already set from X edges)
x_line = np.array([X[0], X[-1]])
# R = (1/thr_best) / T  (in linear space)
y_line = (1.0 / thr_best) / x_line

# Clip to plotted y-range so it doesn't look "out of bounds"
ymin, ymax = Y[0], Y[-1]
mask = (y_line >= ymin) & (y_line <= ymax)
ax1.plot(x_line[mask], y_line[mask], linewidth=3, linestyle="--", color="m",
         label=rf"Threshold line ($\eta_{{best}}$)")

# Panel B
ax2 = plt.subplot(1, 2, 2)
SI_norm = SI / np.nansum(SI) if np.nansum(SI) > 0 else SI
ax2.plot(K, SI_norm, linewidth=1.5, marker="o", color="b", label="SI")

K_pos = K[NN > 0]
NN_pos = NN[NN > 0]
NN_norm = NN_pos / np.sum(NN_pos) if NN_pos.size > 0 else NN_pos
ax2.bar(K_pos, NN_norm, width=1, align="edge",
        facecolor="green", edgecolor="green", linewidth=1, alpha=0.5, label="H(K)")

ax2.axvline(x=K_best, linewidth=3, color="m", linestyle="--", label=rf"$K_{{best}}={K_best}$")
ax2.set_title("NN Histogram", fontsize=14)
ax2.set_xlabel("K", fontsize=14)
ax2.legend(fontsize=14)

plt.savefig("Fig_1_SI.png", bbox_inches="tight", dpi=400)
