import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


SECONDS_PER_YEAR = 365.25 * 24 * 3600


# -------------------- Load output catalog --------------------
df = pd.read_csv("../results/output_catalog.txt", sep=r"\s+", header=None, engine="python")
df.columns = ["time", "lat", "lon", "mag", "bg_id"]
df = df.astype(float)

bg = df[df["bg_id"] == 0]
main = df[df["mag"] >= 6]

t_all_years = df["time"].values / SECONDS_PER_YEAR
t_bg_years = bg["time"].values / SECONDS_PER_YEAR

cum_all = np.arange(1, len(df) + 1)
cum_bg = np.arange(1, len(bg) + 1)

prc_bg = np.round((len(bg) / len(df)) * 100, 2) if len(df) > 0 else 0.0


# -------------------- Plot --------------------
plt.rcParams.update({"font.size": 14})

fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True, dpi=400)
fig.patch.set_facecolor("white")

# Panel A: initial catalog
ax1 = axes[0]
ax1.set_facecolor((0.7, 0.7, 0.7))
ax1.plot(t_all_years, df["lat"].values, "o", markersize=2, markerfacecolor="cyan", markeredgecolor="none")
ax1.plot(main["time"].values / SECONDS_PER_YEAR, main["lat"].values, "o",
         markersize=8, markerfacecolor="red", markeredgecolor="black")

ax1.set_ylabel("Latitude")
ax1.set_title("Initial Catalog")

ax1r = ax1.twinx()
ax1r.plot(t_all_years, cum_all, linewidth=1.5, color="m")
ax1r.set_ylabel("Cum. Num.")

# Panel B: background (declustered)
ax2 = axes[1]
ax2.set_facecolor((0.7, 0.7, 0.7))
ax2.plot(t_bg_years, bg["lat"].values, "o", markersize=1, color="k")
ax2.plot(main["time"].values / SECONDS_PER_YEAR, main["lat"].values, "o",
         markersize=8, markerfacecolor="red", markeredgecolor="black")

ax2.set_title(f"Susceptibility-Based Declustering ({prc_bg}% background events)")
ax2.set_xlabel("Time since origin (years)")
ax2.set_ylabel("Latitude")

ax2r = ax2.twinx()
ax2r.plot(t_bg_years, cum_bg, linewidth=1.5, color="m")
ax2r.set_ylabel("Cum. Num.")

plt.savefig("Fig_2.png", bbox_inches="tight", dpi=400)
