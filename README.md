# Follicle Rupture

Small MATLAB shell simulation for actomyosin-driven follicle rupture.

Run:

```matlab
run_simulation
```

The simulation writes `results/follicleShellSim_test2.mat`. Set
`runPostprocess = true` in `src/follicle_rupture_msd_fixed.m` to generate
diagnostic plots immediately after the run.
Set `makeRuptureMotionVideo = true` to also export
`results/follicleRupture_motion.mp4`.
The movie renderer lives in `src/export_follicle_rupture_motion_video.m` and
can also be called directly on a saved `.mat` file.

Source lives in `src/`.

The passive cell-cell connection model is selected in
`src/follicle_rupture_msd_fixed.m` with:

```matlab
param.connectionModel = 'burgers';
```

Available values are `spring`, `damper`, `maxwell`, `kelvin-voigt`, `zener`,
and `burgers`. Per-edge models can be supplied with `param.connectionModel_edge`
as a cell or string array.

Low-power runtime knobs live near the top of `src/follicle_rupture_msd_fixed.m`:

```matlab
numSave = 300;
runPostprocess = false;
makeRuptureMotionVideo = false;
odeSolver = 'ode45';
progressEveryPhysSec = 5;
```

Try `odeSolver = 'ode15s'` or `'ode23t'` if the damped viscoelastic system
appears stiff on your machine.
