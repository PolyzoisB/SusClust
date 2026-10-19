# Calendar catalog converter

`cat10_to_cat4.f90` converts calendar timestamps to elapsed seconds and removes the depth column. 

## Inputs and outputs

Run from a working directory containing `datasets/` and `results/`. 

- Input: `datasets/input_catalog.txt`.
- Output: `results/output_catalog.txt`, replaced on every successful conversion.

Input must be chronological, headerless, and whitespace-separated:

| Column | Value |
| --- | --- |
| 1-5 | Integer year, month, day, hour, minute |
| 6 | Seconds, including fractional seconds |
| 7 | Signed longitude in degrees |
| 8 | Latitude in degrees |
| 9 | Depth (read, then discarded) |
| 10 | Magnitude |

```text
2021 1 1 0 0 1.25 -118.29 36.05 8.0 3.13
2021 1 1 0 1 2.50 -118.30 36.06 9.0 3.20
2021 1 1 0 2 3.75 -118.31 36.07 7.0 3.40
```

Output columns are **elapsed_seconds latitude longitude magnitude**, matching the SI reader's coordinate order. 

Elapsed time is measured from the first input event, so the first output time is zero. All rows are retained, up to the configured `max_events` limit (2,000,000 by default); empty input stops with an error before writing output. The program does not sort events or validate calendar dates. Use valid Gregorian dates (years at least 1), one consistent time zone, and ordinary seconds; time zones and leap seconds are not handled.

## Build and run

From the repository root in Bash:

```bash
cd Preprocess_tools
mkdir -p datasets results
# Place your ten-column catalog at datasets/input_catalog.txt.
gfortran -O3 cat10_to_cat4.f90 -o cat10_to_cat4.out
./cat10_to_cat4.out
```

In PowerShell:

```powershell
Set-Location Preprocess_tools
New-Item -ItemType Directory -Force -Path datasets,results | Out-Null
# Place your ten-column catalog at datasets/input_catalog.txt.
gfortran -O3 cat10_to_cat4.f90 -o cat10_to_cat4.exe
.\cat10_to_cat4.exe
```

GNU Fortran 10.3 successfully built this program. The converter has no command-line arguments. Edit the source and rebuild to change paths or capacity. Keep its working directory separate from the SI tool's directories: both programs use the same relative input and output filenames for different schemas.

## Use the output with SI

The four-column catalog supplied in the SI `datasets/` directory is already converted and does not need this program. For a new calendar catalog, copy the converter output to `Susceptibility Index/Susceptibility Index tool/datasets/input_catalog.txt`, then follow the [declustering guide](../Susceptibility%20Index/Susceptibility%20Index%20tool/README.md).