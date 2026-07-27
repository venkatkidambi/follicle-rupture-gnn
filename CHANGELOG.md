# Changelog

All notable changes to this project will be documented in this file.

## v0.2.1-bare-bones - 2026-07-27

### Added

- Added centralized cell-scale coordinate configuration with selectable
  `cells`, `microns`, and `meters` input modes.
- Added selectable `cells` or `microns` display units for mesh logging,
  diagnostic plots, and motion videos.
- Added complete SI and display metadata in `param.units`.
- Added lightweight MATLAB unit-conversion validation.

### Changed

- Expressed the default generated ellipsoid as 50 × 15 × 15 cell-diameter
  semi-axes while preserving its 500 µm × 150 µm × 150 µm physical geometry.
- Converted generated and experimental coordinates to meters once at geometry
  ingestion.
- Saved `L0_edge0` as a top-level MAT variable in addition to its existing
  location in `param`, with all physical coordinates, displacements, lengths,
  and times remaining in SI units.
- Centralized plotting conversions instead of repeating hard-coded micron
  scale factors.
- Documented coordinate inputs and the SI units required by every mechanical
  parameter category.

### Compatibility

- Older MAT files without `param.units` continue to be interpreted internally
  in meters and displayed in microns.

### Validation

- Verified 10 µm and 12 µm cell-diameter conversions, equivalent micron/meter
  inputs, cell/micron display conversions, invalid configurations, and legacy
  metadata fallback in MATLAB.

### Known Concerns

- The nondimensional acceleration wrapper may algebraically cancel `mVal`.
  This release reports the concern but does not change the mechanics or
  nondimensionalization equations.

## v0.2.0-bare-bones - 2026-06-29

### Changed

- Reworked the simulation into a bare-bones passive rupture baseline.
- Kept the ellipsoidal epithelial shell, passive viscoelastic edge network,
  localized rupture damage, edge weakening, and edge failure workflow.
- Updated dashboard postprocessing to focus on rupture damage and shell motion
  rather than actomyosin recruitment diagnostics.
- Updated the motion-video exporter for the stripped simulation state.
- Rewrote the README around the current baseline model, quick-start workflow,
  rupture parameters, and future shortening-extension notes.

### Removed

- Removed active actomyosin force generation.
- Removed active rest-length shortening.
- Removed adaptive neighbor-based recruitment and shortening spread.
- Removed the `A_edge` recruitment state from the ODE state.

### Notes

- This version is intentionally minimal: it provides a clean
  rupture-and-weakening baseline that can be extended later with a new
  shortening mechanism.

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
