# Susceptibility Index for Earthquake Declustering

This program computes the Susceptibility Index for an earthquake catalog using the Baiesi–Paczuski metric
and a smoothing-based minimum detection algorithm.

---
## Project Structure

├── global_parameters.f90   # Global free parameters
├── main_program.f90        # Main program
├── randgen.f               # Random number generator
├── datasets/
│ └── input_catalog.txt     # Input catalog
├── results/                # Output files (auto-generated)

---

## Requirements

- GNU Fortran (gfortran)
- Tested with: gfortran 9.4.0

---

## Compilation
From the project root directory, run:
gfortran -O3 global_parameters.f90 main_program.f90 randgen.f -o susceptibility.out

## Execution
./susceptibility.out

## Input Format
The input file must be located at:
datasets/input_catalog.txt

Each line must contain:
time(seconds), longitude(deg), latitude(deg), magnitude

## Output Files
Generated in the results/ folder:
1. rescaled_distances.txt --> maxdt(i), maxdr(i)
These are the Baiesi–Paczuski rescaled distances associated with the strongest link for each event.
This file contains one row per event (except first).
Meaning:
maxdt(i) → Rescaled temporal distance
maxdr(i) → Rescaled spatial distance
Units:
maxdt	years × 10^(-bM/2)
maxdr	km^df × 10^(-bM/2)
Purpose:
Used for diagnostic analysis of clustering structure.

2. susc_index.txt --> Clusters,  SI,  NN,  Threshold,  ip,  nlinks
Meaning:
Clusters → Number of clusters (kr(ip))
SI →	Susceptibility index
NN → Number of events whose strongest similarity lies in this bin
Threshold	→ Similarity threshold value
ip → Log-bin index
nlinks → Total number of links at that threshold
Units:
Threshold → 1 / (seconds × meters^df × 10^(-bM))
Purpose:
Defines the susceptibility curve versus similarity threshold
3. final_results.txt --> kmin_shft, zmin_shft, thresh(kmin_shft), kmin, z_min, thresh(kmin), kmax, thresh(kmax), kmax2, thresh(kmax2)
Meaning:
kmin_shft	→ Number of clusters at shifted minimum
zmin_shft	→ Original SI value at shifted minimum
thresh(kmin_shft) →	Similarity threshold at shifted minimum
kmin	→ Clusters at absolute minimum
z_min	→ Original SI at minimum
kmax	→ Clusters at left peak
kmax2	→ Clusters at right peak

Purpose:
This determines the optimal threshold used for classification
4. output_catalog.txt --> time (seconds),  longitude,  latitude,  magnitude,  id
Units:
id → 0 = background, 1 = triggered
**This is the final classified catalog**
