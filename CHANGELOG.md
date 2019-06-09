# Changelog

## Unreleased

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
