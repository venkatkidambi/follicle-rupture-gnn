# Saved tuning results

## Fixed-surface position and velocity run

The full joint-objective run is in
[fixed_position_velocity_tuning/](fixed_position_velocity_tuning/). It starts
from the approved [fixed-surface resistance fit](archive/02_fixed_surface_resistance/)
and retains its static geometry. Position and interval-velocity losses have
equal weights after normalization by training-only observed-motion RMS scales.

- [Position/velocity comparison table](fixed_position_velocity_tuning/position_velocity_comparison.csv)
- [Velocity comparison, training and held-out cells](fixed_position_velocity_tuning/velocity_comparison.png)
- [Trajectory comparison](fixed_position_velocity_tuning/trajectory_comparison.png)
- [Network and rupture-patch comparison](fixed_position_velocity_tuning/network_comparison.png)
- [Report and numerical verification](fixed_position_velocity_tuning/tuning_report.txt)
- [Simulation video](follicle_rupture_motion_fixed_position_velocity_tuned.mp4)
- [Saved fit](fixed_position_velocity_tuning/best_fit.mat)
- [Forward simulation](follicleShellSim_fixed_position_velocity_tuned.mat)
- [Search log](fixed_position_velocity_tuning/optimization.log)
- [Code snapshot](fixed_position_velocity_tuning/source_snapshot.zip) and
  [run provenance](fixed_position_velocity_tuning/run_provenance.json)

All 240 coarse trials completed, followed by six fine-solver candidates
including the source. Coverage: 36 graph proposals, 59 patch-center proposals,
13 patch-size proposals, 48 scalar proposals covering all ten controls, and
83 local edge swaps. Both the selected graph and weak-link pairs changed.
The joint training objective fell from 1.337129 to 1.273694 (4.74%).
Position RMSE changed from 14.067942 to 13.752792 µm on training samples and
16.747194 to 16.471698 µm on held-out samples. Velocity RMSE changed from
0.0736073 to 0.0717694 µm/s on training intervals and 0.0836529 to 0.0821210
µm/s on held-out intervals. The predicted rupture time is 839.900947 s.

This remains a bounded within-specimen fit. Interval velocities use only
endpoints within each split (about 81.8 s apart), and are derived from the
same positions rather than independent measurements. Late body-band motion
remains substantially underestimated. Tighter integration changed node
positions by at most 0.0000512 µm and rupture time by 0.00000804 s. The
221-frame video was decoded and verified through the final simulation time.

## Moving-surface results

Start with [current/](current/) for the latest full moving-surface refit.
Each run keeps its fit, parameters, search history, reports, simulation,
plots, and video together. Earlier stages are preserved in [archive/](archive/),
and experimental reference videos are in [experimental/](experimental/).

```text
results/
  README.md
  run_summary.csv
  organization_manifest.csv
  current/                              latest full moving-surface refit
  archive/
    01_fixed_surface_network/           original network/mechanics fit
    02_fixed_surface_resistance/        resistance refinement
    03_moving_surface/                  moving surface, fixed-network refit
  experimental/                        tracked experimental reference videos
```

## Latest result

- [Simulation video without pink rupture highlights](current/follicle_rupture_motion_moving_network_tuned_no_rupture_highlights.mp4)
- [Original video with rupture highlights](current/follicle_rupture_motion_moving_network_tuned.mp4)
- [Comparison with previous fits](current/moving_surface_comparison.png)
- [Network and rupture-patch comparison](current/network_comparison.png)
- [Node snapshots](current/node_position_snapshots_moving_network_tuned.png)
- [Ring velocity](current/percentile_ring_velocity_moving_network_tuned.png)
- [Rupture and retraction forces](current/rupture_retraction_force_moving_network_tuned.png)
- [Full report](current/tuning_report.txt)
- [Parameters](current/parameters.csv) and [connections](current/connections.csv)
- [Saved fit](current/best_fit.mat) and [forward simulation](current/follicleShellSim_moving_network_tuned.mat)
- [Test report](current/test_report.txt) and [numerical verification](current/verification.mat)
- [Search checkpoint](current/search_history.mat) and [search log](current/moving_network_search.log)

The video renderer accepts `struct('showRuptureHighlights', false)` as its
third argument to hide the pink rupture edges, circles, and crosses. The
default is `true`. Both videos use the same saved simulation and playback
timing; the new render's metadata is in
[video_info_no_rupture_highlights.mat](current/video_info_no_rupture_highlights.mat).

## Retained tuning stages

The numbers below are copied from each run's saved report. The
[machine-readable summary](run_summary.csv) also lists simulation, video,
and report paths and the predicted rupture time in seconds.

| Stage | What was fitted | Training RMSE (µm) | Reported held-out RMSE (µm) | Report |
| --- | --- | ---: | ---: | --- |
| [01 — Fixed-surface network](archive/01_fixed_surface_network/) | Connections, rupture patch, and seven global controls | 14.09 | 16.75 | [Report](archive/01_fixed_surface_network/tuning_report.txt) |
| [02 — Fixed-surface resistance](archive/02_fixed_surface_resistance/) | Three resistance controls, with the previous network retained | 14.07 | 16.75 | [Report](archive/02_fixed_surface_resistance/tuning_report.txt) |
| [03 — Moving surface](archive/03_moving_surface/) | Ten global controls, with the previous network and patch retained | 13.41 | 12.42 | [Report](archive/03_moving_surface/tuning_report.txt) |
| [04 — Full moving-surface refit](current/) | Network, rupture-patch location/size, and all ten global controls | 12.70 | 11.94 | [Report](current/tuning_report.txt) |

Moving-surface scores are conditional on measured geometry at all frames.
These are comparisons within one specimen, not independent temporal forecasts.
The original tracked [parameter preset](../tuned_parameters.csv) remains at
the repository root because existing scripts read it there.

## Experimental reference

- [Raw tracked motion](experimental/experimental_tracked_motion_raw.mp4)
  and [preview](experimental/experimental_tracked_motion_raw_preview.png)
- [Pose-corrected tracked motion](experimental/experimental_tracked_motion_registered.mp4)
  and [preview](experimental/experimental_tracked_motion_registered_preview.png)
- [Reference metadata](experimental/experimental_reference_info.mat)

## Paths and preservation

The old tuning-directory names and loose result filenames in `results/`
are relative **compatibility symlinks** to the folders above. They are not
duplicate results. They preserve paths used by the unchanged MATLAB scripts,
earlier documentation, and saved metadata. Browse `current/`, `archive/`, and
`experimental/` rather than the compatibility aliases.

The [organization manifest](organization_manifest.csv) records every original
file's old path, consolidated path, byte count, and SHA-256 digest. All 72
original files, including the two root housekeeping files, were preserved and
verified against that inventory. Original result contents and names were not
rewritten. The 44 MATLAB files and the tracked parameter preset were also
verified unchanged during consolidation.

The result artifacts and local compatibility symlinks remain ignored by Git.
This index, the summary table, and the preservation manifest are not ignored
and can be committed as lightweight metadata. A fresh clone does not include
the local simulation artifacts.
