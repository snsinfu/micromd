# Changelog

## Unreleased

### New features

- Forcefield templates:
  - Added `plane_surface_forcefield`: Forcefield on both sides of a plane.
- Potentials:
  - Added `semispring_potential`: Spring potential that is only effective for
    distances longer than the equilibrium distance.
  - Added `wrapped_potential`: Identity wrapper for potential functor. This
    introduces custom potential in the `md` namespace, which enables arithmetic
    operations on the custom potential.

### Changes

- Renamed `system::require()` to `system::add_attribute()` for API consistency.
  The old function `system::require()` is kept for backward compatibility for
  now; it will be deprecated in the future.


## v0.4.0

### Optimizations

- `neighbor_pair_forcefield_v2` is now ~20% faster in some case

### New features

- Potentials:
  - Binary potential functors now support linear arithmetics

## v0.3.1

## v0.3.0

### New features

- Forcefield templates:
  - `neighbor_pair_forcefield_v2`
- Misc:
  - `open_box`, `periodic_box`, `xy_periodic_box`
  - `neighbor_searcher_v2<Box>`

## v0.2.2

### Bug fixes

- Worked around false-positive stack-buffer-overflow reports
- `system::view` performance hit due to unnecessary virtual calls

## v0.2.1

### Bug fixes

- `system::view(key)` is now almost zero-overhead

## v0.2.0

### New features

- Forcefield templates:
  - `intra_subsystem_pair_forcefield`
  - `inter_subsystem_pair_forcefield`
  - `intra_subsystem_neighbor_pair_forcefield`
  - `inter_subsystem_neighbor_pair_forcefield`

### Changes

- Renamed `subsystem_pair_forcefield` to `intra_subsystem_pair_forcefield`

## v0.1.1

### New features

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

### New features

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

### New features

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
