# ETAS simulations and SI evaluation

This program is intended to sample ETAS parameter sets, generate synthetic catalogs, apply Susceptibility Index (SI) declustering, and compare the result with known background/triggered labels.

## Build blocker in this checkout

**The required `global_parameters.f90` is missing from this directory.** `main_program.f90` imports that module, and the former README's build command references it. The module in the neighboring SI tool directory does not provide the ETAS-specific variables and is not a substitute. Restore the ETAS configuration module before building or running this program.

This guide describes the behavior present in the source, not a verified runnable ETAS workflow. Default simulation counts, parameter ranges, windows, paths, and seeds cannot be stated reliably without the missing module.

## Files present

```text
Susceptibility Index/ETAS simulations/
  main_program.f90
  randgen.f
  datasets/
    cat_main_m2.dat
    cat_main_m4.dat
  results/
    combinations.txt
    ssm_min_final.txt
```

The two checked-in result files contain no simulation records. Detailed diagnostics are written only for simulations outside the configured ratio band; those directories are not included in this checkout.

## Simulation workflow

For each of `num_sets` simulations, the source:

1. Samples parameters within configured ranges, rejecting combinations outside the branching-ratio interval.
2. Samples background coordinates from two magnitude-dependent pools, adds spatial jitter, and generates parent events and their offspring in an auxiliary region.
3. Sorts events by time, tracks parent relationships, and selects a target time/space/magnitude window.
4. Computes SI over similarity thresholds and selects an unshifted and shifted minimum.
5. Compares estimates with the number of true background events retained in the target window and evaluates background-event classification.
6. Writes parameters and summary metrics; exports detailed curves and a classified target catalog if `K_best / y_true` is outside the configured band.

`y_true` counts target events with `parent_id = 0`. It is not a count of every distinct ancestral cluster represented in the target window: descendants can have roots outside that window.

## Background coordinate pools

`load_bg_coords()` reads six numeric, whitespace-separated columns from each pool, with no header:

| Column | Use in the source |
| --- | --- |
| 1 | Read and ignored |
| 2 | Read and ignored |
| 3 | Magnitude; rows from the lower-magnitude pool are included only if magnitude is at least 1 |
| 4 | `x`, used as longitude in the simulation and metric |
| 5 | `y`, used as latitude in the simulation and metric |
| 6 | Read and ignored |

Locations for sampled parent magnitudes at least 4 are drawn from the `m4` pool; other parent locations come from the retained `m2` pool. Both pools are already supplied in `datasets/`; their `.dat` files are identical to the corresponding `.txt` files in `Data/Simulations/Southern California/Background coordinates/`.

**Reconcile the coordinate convention before using them.** The bundled fourth column ranges roughly from 30 to 37.5 and the fifth from 113 to 122, consistent with latitude followed by positive west longitude for Southern California. The current ETAS code assigns those columns to longitude and latitude respectively. After confirming that convention against the data source, put signed longitude in column 4 and latitude in column 5 (for these western coordinates, longitude should be negative), or adapt the loader and all downstream coordinate handling consistently. Match the restored spatial-window parameters to the same convention. Merely changing file extensions does not resolve this mismatch.

## Configuration to restore

The missing module must define the following names used by `main_program.f90`. Recover the intended scientific values rather than copying the SI tool's defaults.

| Group | Names referenced by the source |
| --- | --- |
| Capacity | `max_bg_in`, `max_events`, `q_buf_size`, `max_thr`, `max_points` |
| Files/directories | `bg_coords_m4_file`, `bg_coords_m2_file`, `results_dir`, `ssm_index_dir`, `sim_catalog_dir`, `summary_file`, `combinations_file` |
| Simulation count and seed | `num_sets`, `seed_default` |
| Time conversion and window | `seconds_per_day`, `seconds_per_year`, `tc_years`, `tlast_years` |
| Auxiliary spatial window | `x0min`, `x0max`, `y0min`, `y0max`, `coord_jitter_deg` |
| Target spatial/magnitude window | `xmin`, `xmax`, `ymin`, `ymax`, `qc` |
| Simulated magnitude range | `qmin`, `qsup`, `beta_val` |
| Parameter ranges | `pp_min/max`, `cc_min_days/max_days`, `alpha_min/max`, `kk_min/max`, `gamma_min/max`, `dd_min/max`, `qdec_min/max`, `bg_rate_min/max`, `nbranch_min/max` |
| SI scan and metric | `pthmin`, `pthmax`, `bin0`, `smooth_step`, `df`, `bval` |
| Diagnostic export band | `export_ratio_lo`, `export_ratio_hi` |

Names written as `pp_min/max` mean the two identifiers `pp_min` and `pp_max`; similarly, `cc_min_days/max_days` means `cc_min_days` and `cc_max_days`.

Time inputs `tc_years` and `tlast_years` are converted to seconds. The output catalog stores seconds from simulation origin. `cc_days` is recorded in days before conversion for simulation; the configured `dd` values and coordinate windows must match the spatial convention of the simulation. SI similarities use `1/(seconds * km**df)` here, whereas the standalone SI tool uses years. Rescaled `T` uses years and `R` uses `km**df`.

## Build and run after restoring the module

Use GNU Fortran and a POSIX shell such as Bash. From the repository root, **after restoring `global_parameters.f90`, reconciling the pools, and configuring the paths**:

```bash
cd "Susceptibility Index/ETAS simulations"
gfortran -O3 global_parameters.f90 main_program.f90 randgen.f -o etas_ssm.out
./etas_ssm.out
```

Run from this directory so relative paths resolve. The source calls `mkdir -p` to create output directories; this assumes a POSIX shell and is not portable to native Windows. For Windows, use WSL or adapt `ensure_dir()` before attempting a native run. Summary and parameter files are replaced on each run; save previous results separately.

## Output schemas

The actual output paths come from the restored module. The following names follow the existing result files and the source's diagnostic naming convention. Files contain numeric whitespace-separated rows, without headers.

### Parameter combinations

One row per simulation, in this exact order:

```text
alpha kk pp cc_days DD gamma qdec bg_rate beta n_branch
```

### Summary evaluation

`summary_file` receives 17 columns:

| Column | Value |
| --- | --- |
| 1 | `n_c`, retained target events |
| 2 | `y_true`, true target background events |
| 3 | `f1_true`, F1 at the curve position where `K` crosses `y_true` (0 if no crossing is found) |
| 4 | `k_max`, cluster estimate at maximum F1 |
| 5 | `f1_max`, maximum F1 over scanned curve points |
| 6 | `kmin`, unshifted-minimum cluster estimate |
| 7 | `kmin_shft` (`K_best`), shifted-minimum cluster estimate |
| 8 | `f1_min`, F1 at the unshifted minimum |
| 9 | `f1_shft`, F1 at the shifted minimum |
| 10 | `prec_shft`, background precision at the shifted minimum |
| 11 | `rec_shft`, background recall at the shifted minimum |
| 12 | `acc_shft`, classification accuracy at the shifted minimum |
| 13 | True triggered fraction, `(n_c-y_true)/n_c` |
| 14 | Estimated triggered fraction, `(n_c-kmin_shft)/n_c` |
| 15 | `kmin_shft/y_true` |
| 16 | `kmin/y_true` |
| 17 | `k_max/y_true` |

The code's confusion counters treat triggered events as positive internally, then compute precision, recall, and F1 for **background** classification using the true-negative counter. Scores cover target events 2 through `n_c`; fractions and count ratios use all target events. Account for this distinction when interpreting the summaries.

### Conditional detailed diagnostics

If `kmin_shft/y_true < export_ratio_lo` or `> export_ratio_hi`, the program writes:

- `ssm_index_comb_X.txt` under `ssm_index_dir`: 15 columns in the order below.
- `sim_cat_comb_X.txt` under `sim_catalog_dir`: 8 columns in the order below.

```text
# SI diagnostic columns (no header is actually written):
K SI nd nd_bg nd_tr nd/n_c nlink/n_c derivative Threshold ThresholdIndex Nlinks precision recall F1 accuracy

# Target catalog columns (no header is actually written):
time_seconds longitude latitude magnitude event_id parent_id truth_cluster_id predicted_label
```

`parent_id = 0` identifies a true background event. `predicted_label = 0` means predicted background, and `1` means predicted triggered; the first target event is always labeled 0. ETAS catalog coordinates are longitude before latitude, unlike the standalone SI tool's catalog. `X` is the simulation iteration (for example, `comb_1`). The parameter and summary rows correspond by iteration, as neither file includes an explicit iteration column.

## Further implementation limits

- A single configured seed does not control all randomness: event sampling uses `ran2`, parameter sampling uses Fortran `random_number`, and Poisson sampling uses `randgen.f`. The program does not explicitly seed the latter two streams. Record and initialize each stream before expecting repeatable simulations.
- `cluster_index()` stops at `N-1`, leaving the final event's truth cluster ID unset. Repair this before relying on exported cluster IDs.
- Event and pool arrays have fixed capacities without comprehensive overflow checks. Pairwise work and threshold-bin updates can be expensive; link counters can overflow.
- The warning for similarities outside the threshold range does not prevent invalid array indexing. Empty target catalogs, zero true background counts, and SI curves without suitable extrema lack robust handling.

These issues remain source-code work. See the [documentation review](../../DOCUMENTATION_REVIEW.md) and [root citation section](../../README.md#citation).
