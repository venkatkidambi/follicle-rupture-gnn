# Changelog

All notable changes to this project will be documented in this file.

## v0.1.0-alpha - 2026-06-29

### Added

- Initial alpha release of the follicle rupture shell simulation.
- MATLAB wrapper entry point in `run_simulation.m`.
- Main ellipsoidal epithelial shell model in `src/follicle_rupture_msd_fixed.m`.
- Active rest-length shortening mode for actomyosin-like contraction.
- Adaptive active recruitment dynamics with posterior tip seeding.
- Optional active tension model for comparison.
- Dashboard postprocessing in `src/confirm_actomyosin_working.m`.
- MP4 motion export in `src/export_follicle_rupture_motion_video.m`.
- Experimental node data loading workflow for spreadsheet-based geometry tests.
- Example input-node spreadsheet and generated result artifacts.

### Notes

- This is an experimental prerelease, not a stable full release.
- Experimental node data support remains work in progress and should be treated
  as diagnostic until geometry, edge graph, force projection, and parameter
  scales are validated for each dataset.
