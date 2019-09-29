# `THDMMinimizer.jl` Scripts

This directory contains various scripts used to find simultaneous charge-breaking and electo-weak-like minima at 1-loop in the Two-Higgs-Doublet (2HDM) model. Below, we give brief descriptions of the various scripts.

## `Scanner.jl`
This is the main file for finding 2HDM parameters which yield simultaneous charge-breaking and electro-weak minima. The algorithm used to discover these parameters is laid out in detail in the paper. A brief overview of the algorithm is:

1. generate a normal and CB vacuum and THDM parameters such that the tree-level potential is approximately extremized at the vacuua and the tree-level potential is bounded.
2. Choose 50 new randomly placed vacuua and perform a minimization at each in order to possibly find all minima.
3. Determine if parameters yeild a normal and CB minimum.
4. Determine which is deepest
5. Save data.

The algorithm can be ran by first installing the package or making the package visible to `julia`, navigating to the `scripts` (this) directory, starting a `julia` session and running:
```julia
julia> include("scanner.jl")
# Discover only points which have simultaneous CB and EW minima
julia> scan(onlya=true)
# Discover points which have simultaneous CB and/or EW minima
julia> scan()
```
A colored promt will appear and display the current number of points that have been found. The results will be saved in the directory in which you ran the script. For points where there are CB and normal minima, the results will be saved in a file called `typea1.csv` (if the CB minima is deepest) and `typea2.csv` (if the normal is deepest). You one set `onlya=true` or omitted the keyword argument, points in which there are only EW minima will be stored in a file named `typeb.csv` and if there is only a CB minima, stored in `typec.csv`. This script will run endlessly and will need to be stopped by the user.

## `scan_from_nuclei.jl`
This script uses the values from `data/verified_a1.csv` and `data/verified_a2.csv` as starting points, branches off to nearby vacuua and tries to find new counterexamples. It can be ran following the same steps as above and running:
```julia
# for type A1 points
julia> scan(:a1)
# for type A2 points
julia> scan(:a2)
```

## `data_verification.jl`
This script will run through all of the points in `data/candidates_a1.csv`, `data/candidates_a2.csv`, etc. and determine the points which are perturbative, bounded at tree-level and have the correct number of Golstone bosons and save them to `verified_a1.csv`, `verified_a2.csv`, etc.

## Plotting

### `1Dplots.jl`
This script will generate the 1D-plots used in the paper. It can be modified to create plots for any of the points one wishes to plot.

### `2Dplots.jl`
This script will generate the 2D-plots used in the paper. It can be modified to create plots for any of the points one wishes to plot.

### `rge_plots.jl`
This script will generate the RG plots used in the paper.
It performs the running of the parameters and minimizations for various renormalization scales. It can be modified to create plots for any of the points one wishes to plot.
