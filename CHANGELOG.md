# Changelog

## Unreleased

### Added

- Forcefield templates:
  - `neighbor_pair_forcefield_v2`
- Misc:
  - `open_box`, `periodic_box`, `xy_periodic_box`
  - `neighbor_searcher_v2<Box>`


## v0.2.2

### Fixed

- Worked around false-positive stack-buffer-overflow reports
- `system::view` performance hit due to unnecessary virtual calls

## v0.2.1

### Fixed

- `system::view(key)` is now almost zero-overhead

## v0.2.0

### Added

- Forcefield templates:
  - `intra_subsystem_pair_forcefield`
  - `inter_subsystem_pair_forcefield`
  - `intra_subsystem_neighbor_pair_forcefield`
  - `inter_subsystem_neighbor_pair_forcefield`

### Changed

- Renamed `subsystem_pair_forcefield` to `intra_subsystem_pair_forcefield`

## v0.1.1

### Added

- Forcefield templates:
  - `selected_pair_forcefield`
  - `subsystem_pair_forcefield`
- Potentials:
  - `polybell_potential`
- Misc:
  - `particle_ref::view(key)` for attribute access
  - `power<N>(x)`
  - `power_sqrt<N>(x)`

### Fixed

- `sphere_surface_forcefield` overwrites force array

## v0.1.0

### Added

- Forcefield templates:
  - `all_pair_forcefield`
  - `point_source_forcefield`
  - `sequential_triple_forcefield`
- Potentials:
  - `soft_lennard_jones_potential`
  - `spring_potential`
- Misc:
  - `neighbor_searcher`

## v0.0.1

### Added

- Basic functionalities
- Forcefield templates:
  - `neighbor_pair_forcefield`
  - `sequential_pair_forcefield`
  - `sphere_surface_forcefield`
  - `ellipsoid_surface_forcefield`
- Potentials:
  - `constant_potential`
  - `harmonic_potential`
  - `lennard_jones_potential`
  - `softcore_potential`
- Simulations:
  - `brownian_dynamics_simulation`
  - `newtonian_dynamics_simulation`
