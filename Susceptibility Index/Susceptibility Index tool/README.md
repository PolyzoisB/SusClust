# Susceptibility Index tool

This tool reads an earthquake catalog, computes nearest-neighbor similarities and the Susceptibility Index (SI), selects a threshold using a smoothed and shifted SI minimum, and labels retained events as background or triggered. It exports binary event labels, rather than a cluster assignment for every earthquake.

## Build and run

Use GNU Fortran. From the repository root in Bash:

```bash
cd "Susceptibility Index/Susceptibility Index tool"
mkdir -p results
gfortran -O3 global_parameters.f90 main_program.f90 -o susceptibility.out
./susceptibility.out
```

In PowerShell:

```powershell
Set-Location 'Susceptibility Index/Susceptibility Index tool'
New-Item -ItemType Directory -Force -Path results | Out-Null
gfortran -O3 global_parameters.f90 main_program.f90 -o susceptibility.exe
.\susceptibility.exe
```

The build succeeds with GNU Fortran 10.3; the legacy `PAUSE` statement produces a warning. `randgen.f` is included in the command for consistency with the supplied sources, but no random-number routines are called by this tool.

Keep the working directory at `Susceptibility Index/Susceptibility Index tool` when running the executable. `datasets/` and `results/` are relative to that directory. Create `results/` before running: the program does not create it. Each run replaces the four output files. There are no command-line arguments; configure the source and recompile.

## Input catalog

`datasets/input_catalog.txt` must contain four whitespace-separated numeric columns, with no header:

| Column | Meaning | Units/convention |
| --- | --- | --- |
| 1 | Elapsed time from a chosen origin | Seconds |
| 2 | Latitude | Degrees |
| 3 | Longitude | Signed degrees; west is negative |
| 4 | Magnitude | Earthquake magnitude |

```text
140589.219 36.04838 -118.29092 3.13
```

## Configuration

Edit `global_parameters.f90` and rebuild:

| Parameter | Default | Effect |
| --- | --- | --- |
| `input_file` | `datasets/input_catalog.txt` | Input path relative to the working directory |
| `mcut` | `3.0` | Inclusive magnitude cutoff |
| `df` | `1.6` | Spatial exponent in the nearest-neighbor metric |
| `bval` | `1.0` | Fixed Gutenberg-Richter b-value used by the metric |
| `pthmin`, `pthmax` | `1e-18`, `1e15` | Similarity threshold scan limits |
| `bin0` | `0.05` | Natural-log threshold spacing |
| `max_events` | `2000000` | Input-row reading limit and event-array capacity |
| `max_thr` | `6000` | Index capacity of threshold arrays, from `-max_thr` to `max_thr` |
| `max_points` | `10000` | SI curve-array capacity |
| `rescaled_dist`, `susc_index_file_name`, `output_cat`, `final_results` | Paths under `results/` | Output destinations |


## Metric and units

For an event and an earlier event of magnitude `m`, the similarity is

```text
s = 1 / (dt_years * dr_km**df * 10**(-bval*m))
T = dt_years * 10**(-bval*m/2)
R = dr_km**df * 10**(-bval*m/2)
T * R = 1 / s
```

## Output files

All four files are whitespace-separated and headerless when produced by the current source.

| File | Contents |
| --- | --- |
| `results/rescaled_distances.txt` | Two columns, `T R`, for each retained event after the first. The values come from the link that first reaches that event's highest similarity bin; within-bin ties are not resolved using the exact similarity. |
| `results/susc_index.txt` | Six columns: `K SI NN Threshold ThresholdIndex Nlinks`. `K` is the algorithm's cluster-count estimate, `NN` counts strongest-link bins, `Threshold = exp(ThresholdIndex * bin0)`, and `Nlinks` counts links at that threshold. |
| `results/final_results.txt` | One numeric row with the ten columns listed below. |
| `results/output_catalog.txt` | Five columns: `time_seconds latitude longitude magnitude bg_id`, for the retained input events. `bg_id = 0` means background; `bg_id = 1` means triggered. The first retained event is always labeled background. |

Exact `final_results.txt` column order:

| Column | Meaning |
| --- | --- |
| 1 | `K_best`: cluster estimate at the shifted minimum |
| 2 | Original SI value at `K_best` |
| 3 | Similarity threshold at `K_best` |
| 4 | `K_min`: cluster estimate at the unshifted minimum |
| 5 | Original SI value at `K_min` |
| 6 | Similarity threshold at `K_min` |
| 7 | `K_leftPeak` |
| 8 | Similarity threshold at the left peak |
| 9 | `K_rightPeak` |
| 10 | Similarity threshold at the right peak |

The SI values written to this summary come from the original curve at the selected positions, even though smoothing is used to choose those positions. Background classification uses the selected log-bin threshold. Count `bg_id == 0` if you need the actual number of exported background events; the first event's handling means you should not assume that count equals `K_best`.

## Plots

Install NumPy, pandas, Matplotlib, and SciPy using `python -m pip install numpy pandas matplotlib scipy`; see the [root README](../../README.md#requirements) for a virtual environment example. Run from `plots/` after generating fresh Fortran results:

```bash
cd plots
python fig_1.py
python fig_2.py
```

- `fig_1.py` reads `susc_index.txt`, `final_results.txt`, and `rescaled_distances.txt`. It saves `plots/Fig_1_SI.png`, showing a smoothed density of rescaled distances with the selected threshold, and the normalized SI curve and nearest-neighbor histogram.
- `fig_2.py` reads `output_catalog.txt` and saves `plots/Fig_2.png`. It plots latitude against years since origin for retained and background events, with cumulative counts. Red markers show events with magnitude at least 6; that cutoff is in the script, and those markers are overlaid in both panels regardless of background classification.

## Citation

See the [root citation section](../../README.md#citation) and [CITATION.cff](../../CITATION.cff).
