# Post Processing

The scripts in this directory process the model results and produce:

- the tracking files, containing the evolution of the strain rate within the Ocean-Continent Transition (OCT) for all models
- the complete `sr_raw` panels, with every `Nemin` and `SRmean` series for every model, velocity, plasticity and viscosity;
- diagnostic signal/derivative plots directly in `outputs/processing_steps`, named `A_<model>_<velocity>_<viscosity>_<plasticity>.png` or `mSR_<model>_<velocity>_<viscosity>_<plasticity>.png`;
- `tabella_subduzione.csv`, containing the `tSL`, `tSZI`, trends and deltas discussed within the paper. 

Colours, fonts, linewidths and line styles of the figures can be configured in `config.py`.

## Step 1: Strain rate tracking

`SRTracking.f90` extracts the evolution of the strain rate within the OCT from the model results.
The Python workflow (Step 2) expects their `SR2_*.txt` output under `../Tracking/`.

``` bash
    gfortran SRTracking.f90
```

## Step 2: Python workflow

The Python workflow assesses the time of strain rate localisation and subduction zone initiation starting from the tracking files, which are the output of Step 1.

``` bash
    python3 main.py
```

### config.py

Sets directories, colours, variables and parameters for the post processing of the tracking files. Set the analysis you want to perform. 

### main.py

Loads the modules and calls the functions using the parameters, the directories and the options set in config.py. This is the code that has to be run.

### StrainRate.py

Contains deterministic smoothing, split-time detection and the single-series plotter.
