# Susceptibility Index (SI) for Earthquake Declustering

This folder is a self-contained tool that:
1. Reads an earthquake catalog.
2. Computes the Susceptibility Index (SI) using the Baiesi–Paczuski metric.
3. Automatically selects an optimal similarity threshold (shifted minimum criterion).
4. Outputs a declustered catalog (background vs triggered) and diagnostic files.
5. Provides plotting scripts to reproduce key figures.

## Folder layout

```
Susceptibility Index/
├─ global_parameters.f90
├─ main_program.f90
├─ randgen.f
├─ datasets/
│  └─ input_catalog.txt
├─ results/
│  ├─ susc_index.txt
│  ├─ final_results.txt
│  ├─ output_catalog.txt
│  └─ rescaled_distances.txt
└─ plots/
   ├─ fig_1.py
   └─ fig_2.py
```

## Requirements

### Fortran
- `gfortran` (tested with gfortran 9+)

### Python (for plots only)
- Python 3.9+
- `numpy`, `pandas`, `matplotlib`, `scipy`

Install (recommended in a venv):
```bash
pip install numpy pandas matplotlib scipy
```

## Quickstart

Run everything from **this** directory (`Susceptibility Index/`):

### 1) Compile
```bash
gfortran -O3 global_parameters.f90 main_program.f90 randgen.f -o susceptibility.out
```

### 2) Run
```bash
./susceptibility.out
```

This will read the input catalog from:
- `datasets/input_catalog.txt`

and write outputs to:
- `results/`

### 3) Make plots
```bash
cd plots
python fig_1.py
python fig_2.py
```

Plots are saved into the `plots/` folder.

## Input format

The input file must be located at:

- `datasets/input_catalog.txt`

Each line must contain **four whitespace-separated columns**:

1. `time_seconds` — elapsed time since origin (seconds)
2. `latitude_deg`
3. `longitude_deg`
4. `magnitude`

Example:
```
140589.219 36.04838 -118.29092 3.13
```

Notes:
- Events with magnitude `< mcut` are ignored (`mcut` is defined in `global_parameters.f90`).
- Lines with non-positive times are ignored.

## Outputs (written to `results/`)

### 1) `results/rescaled_distances.txt`
Two columns per event (except the first event):
- `maxdt(i)`: rescaled time distance
- `maxdr(i)`: rescaled spatial distance

These are associated with the **strongest link** of each event (largest similarity).

### 2) `results/susc_index.txt`
Columns:
1. `K` — number of clusters
2. `SI` — susceptibility index value
3. `NN` — histogram count of strongest-link bins
4. `Threshold` — similarity threshold value
5. `ThresholdIndex` — log-bin index (`ip`)
6. `Nlinks` — total number of links at that threshold

### 3) `results/final_results.txt`
One line with:
- `K_best` (shifted-minimum clusters), `SI_best`, `Thr(K_best)`,
- `K_min`, `SI_min`, `Thr(K_min)`,
- `K_leftPeak`, `Thr(K_leftPeak)`,
- `K_rightPeak`, `Thr(K_rightPeak)`

(Exact column ordering matches the Fortran write statement.)

### 4) `results/output_catalog.txt`
Final classified catalog with 5 columns:
1. `time_seconds`
2. `latitude_deg`
3. `longitude_deg`
4. `magnitude`
5. `bg_id` — **0 = background**, **1 = triggered**

## Parameters

Edit `global_parameters.f90`:
- `mcut` (magnitude cutoff)
- `df` (fractal dimension)
- `bval` (Gutenberg–Richter b-value, if you want it fixed)
- `pthmin`, `pthmax`, `bin0` (threshold scan range and bin size)
- `max_events` (max catalog size supported by fixed arrays)

## Notes / common pitfalls

- Run the executable from this folder so the relative paths resolve correctly.
- Large catalogs are computationally expensive because the proximity computation is O(N²).
- If you change the input format, update both the input file and the reading logic in `main_program.f90`.

## Citation
If you use this tool in research, please cite the associated SusClust method/publication (add reference here).
