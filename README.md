# Follicle Rupture

Small MATLAB shell simulation for actomyosin-driven follicle rupture.

Run:

```matlab
run_simulation
```

The simulation writes `results/follicleShellSim_test2.mat` and runs diagnostics.

Source lives in `src/`.

The passive cell-cell connection model is selected in
`src/follicle_rupture_msd_fixed.m` with:

```matlab
param.connectionModel = 'burgers';
```

Available values are `spring`, `damper`, `maxwell`, `kelvin-voigt`, `zener`,
and `burgers`. Per-edge models can be supplied with `param.connectionModel_edge`
as a cell or string array.
