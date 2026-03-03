# ETAS Simulations + Susceptibility Index Declustering (Evaluation Tool)

This folder runs synthetic ETAS simulations and then applies the **Susceptibility Index (SI)** declustering algorithm to evaluate performance against the known ETAS parent/offspring structure.

## What it does

For each simulation (repeated `num_sets` times):

1. Randomly samples an ETAS parameter set from predefined ranges.
2. Simulates an ETAS catalog (parents + offspring) in an auxiliary region.
3. Extracts a target spatio-temporal-magnitude window.
4. Computes SI using the Baiesi–Paczuski similarity metric.
5. Finds the SI minimum and shifted-minimum threshold.
6. Compares predicted number of clusters with the true number of background parents.
7. Writes summary results and (optionally) exports detailed diagnostics for “bad” cases.

## Folder layout (recommended)

```
ETAS_simulations/
├─ global_parameters.f90
├─ main_program.f90
├─ randgen.f                 # if you need ZBQLPOI or other RNG utilities
├─ datasets/
│  ├─ cat_main_m2.dat        # background coordinates pool (lower magnitudes)
│  └─ cat_main_m4.dat        # background coordinates pool (higher magnitudes)
└─ results/
   ├─ combinations.txt
   ├─ ssm_min_final.txt
   ├─ ssm_index/             # created automatically
   └─ sim_catalogs/          # created automatically
```

## Requirements

- GNU Fortran (`gfortran` 9+ recommended)

## Input datasets

You must provide two files in `datasets/`:

- `cat_main_m4.dat`
- `cat_main_m2.dat`

They are used as pools of background coordinates for sampling parent-event locations.

The program expects each line to contain at least:
`z t mag x y z`
(i.e., it reads 6 columns and uses the `x` and `y` columns as coordinates.)

If your files have different columns/order, update `load_bg_coords()` accordingly.

## Compile

From inside `ETAS_simulations/`:

```bash
gfortran -O3 global_parameters.f90 main_program.f90 randgen.f -o etas_ssm.out
```

## Run

```bash
./etas_ssm.out
```

Outputs are written to `results/`.

## Outputs

### `results/combinations.txt`
One row per simulation with the sampled ETAS parameters.

### `results/ssm_min_final.txt`
One row per simulation with summary evaluation metrics, including:
- total number of target events `n_c`
- true number of clusters `y_true`
- predicted clusters from `K_min` and `K_best` (shifted min)
- F1 / precision / recall at selected thresholds (as implemented in code)

### Optional diagnostics
If the ratio `K_best / y_true` is outside a tolerance band (see `export_ratio_lo/hi`), the program exports:

- `results/ssm_index/ssm_index_comb_X.txt` (full SI curve + metrics)
- `results/sim_catalogs/sim_cat_comb_X.txt` (target catalog with truth cluster id + predicted bg/trigger label)

## Configuration

Edit `global_parameters.f90` to change:
- number of simulations (`num_sets`)
- parameter ranges (pp, cc, alpha, etc.)
- spatial/temporal/magnitude windows
- SI scan ranges (`pthmin`, `pthmax`, `bin0`)
- export thresholds for diagnostics

## Notes / recommendations

- This code is computationally heavy (pairwise O(N²) similarity for each target catalog).
- Memory usage is high due to fixed-size arrays. If you want to support larger catalogs or reduce memory, consider refactoring to allocatable arrays sized at runtime.
- Make sure you run the executable from this folder, so relative dataset/result paths resolve.
