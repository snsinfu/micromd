# Quick Reference

## Basic types

Simple aliases:

```c++
using index  = std::size_t;
using step   = std::uint64_t;
using scalar = double;
```

Structs:

```c++
struct point {
    scalar x;
    scalar y;
    scalar z;

    scalar distance(point) const;
    scalar squared_distance(point) const;
};

struct vector {
    scalar x;
    scalar y;
    scalar z;

    scalar dot(vector) const;
    vector cross(vector) const;
    vector hadamard(vector) const;
    scalar norm() const;
    scalar squared_norm() const;
    vector normalize() const;
    vector project(vector) const;
};
```


## System

```c++
class system {
    // Particles
    particle_ref   add_particle(data);
    particle_range particles();

    // Forcefield
    shared_ptr     add_forcefield(ff);
    scalar         compute_energy();
    void           compute_force(forces);

    // Attributes
    void           add_attribute(key);
    array_view     view(key);
    array_view     view_masses();
    array_view     view_mobilities();
    array_view     view_positions();
    array_view     view_velocities();
};

struct basic_particle_data {
    scalar mass;
    scalar mobility;
    point  position;
    vector velocity;
};

struct particle_ref {
    scalar& mass;
    scalar& mobility;
    point&  position;
    vector& velocity;

    auto& view(key);
};
```


## Forcefield

Virtual interface:

```c++
class forcefield {
    virtual scalar compute_energy(system);
    virtual void   compute_force(system, forces);
};
```


## Forcefield templates

### Composite

```c++
class composite_forcefield<FF...> {
};
```


### All pairs

CRTP base class:

```c++
class bruteforce_pairwise_forcefield<Derived> {
    Box  unit_cell(system);
    auto bruteforce_targets(system);
    auto bruteforce_pairwise_potential(system, i, j);
};
```

Basic implementation:

```c++
auto make_bruteforce_pairwise_forcefield(pot);
```


### Neighbor pairs

CRTP base class:

```c++
class neighbor_pairwise_forcefield<Derived, Box> {
    Box    unit_cell(system);
    auto   neighbor_targets(system);
    scalar neighbor_distance(system);
    auto   neighbor_pairwise_potential(system, i, j);
};
```

Basic implementation:

```c++
class basic_neighbor_pairwise_forcefield<Derived, Box> {
    this_t set_unit_cell(box);
    this_t set_unit_cell(box_cb);
    this_t set_neighbor_distance(dist);
    this_t set_neighbor_distance(dist_cb);
    this_t set_neighbor_targets(indices);
    this_t set_neighbor_targets(indices_cb);
};

auto make_neighbor_pair_forcefield<Box>(pot);
```


### Bonded pairs

CRTP base class:

```c++
class bonded_pairwise_forcefield<Derived> {
    this_t add_bonded_pair(i, j);
    this_t add_bonded_range(start, end);
    auto   bonded_pairwise_potential(system, i, j);
};
```

Basic implementation:

```c++
auto make_bonded_pairwise_forcefield(pot);
```


### Bonded triples

CRTP base class:

```c++
class bonded_triplewise_forcefield<Derived> {
    this_t add_bonded_triple(i, j, k);
    this_t add_bonded_range(start, end);
    auto   bonded_triplewise_potential(system, i, j, k);
};
```

Basic implementation:

```c++
auto make_bonded_triplewise_forcefield(pot);
```


### Plane surface

CRTP base class:
```c++
class plane_surface_forcefield<Derived> {
    plane plane(system);
    auto  plane_inward_potential(system, i);
    auto  plane_outward_potential(system, i);
};

struct plane {
    vector normal;
    point  reference;
};
```

Basic implementation:

```c++
class basic_plane_surface_forcefield<Derived> {
    this_t set_plane(plane);
    this_t set_plane(plane_cb);
};

auto make_plane_inward_forcefield(pot);
auto make_plane_outward_forcefield(pot);
```


### Sphere surface

CRTP base class:

```c++
class sphere_surface_forcefield<Derived> {
    sphere sphere(system);
    auto   sphere_inward_potential(system, i);
    auto   sphere_outward_potential(system, i);
};

struct sphere {
    scalar radius;
    point  center;
};
```

Basic implementation:

```c++
class basic_sphere_surface_forcefield<Derived> {
    this_t set_sphere(sphere);
    this_t set_sphere(sphere_cb);
};

auto make_sphere_inward_forcefield(pot);
auto make_sphere_outward_forcefield(pot);
```


### Ellipsoid surface

CRTP base class:

```c++
class ellipsoid_surface_forcefield<Derived> {
    ellipsoid ellipsoid(system);
    auto      ellipsoid_inward_potential(system, i);
    auto      ellipsoid_outward_potential(system, i);
};

struct ellipsoid {
    scalar semiaxis_x;
    scalar semiaxis_y;
    scalar semiaxis_z;
    point  center;
};
```

Basic implementation:

```c++
class basic_ellipsoid_surface_forcefield<Derived> {
    this_t set_ellipsoid(ellipsoid);
    this_t set_ellipsoid(ellipsoid_cb);
}

auto make_ellipsoid_inward_forcefield(pot);
auto make_ellipsoid_outward_forcefield(pot);
```


## Potentials

### Pairwise potentials

```c++
// u(r) = e
struct constant_potential {
    scalar energy;
};

// u(r) = K/2 r^2
struct harmonic_potential {
    scalar spring_constant;
};

// u(r) = K/2 (r - b)^2
struct spring_potential {
    scalar spring_constant;
    scalar spring_length;
};

// u(r) = K/2 (r - b)^2   (r > b)
struct semispring_potential {
    scalar spring_constant;
    scalar spring_length;
}

// u(r) = e ((s/r)^12 - 2 (s/r)^6)
struct lennard_jones_potential {
    scalar epsilon;
    scalar sigma;
};

// u(r) = e ((k+1)/(k+(r/s)^6) - 1)^2 - e
struct soft_lennard_jones_potential {
    scalar k;
    scalar epsilon;
    scalar sigma;
};

// u(r) = e (1 - (r/d)^p)^q   (r < d)
struct softcore_potential<P=3, Q=2> {
    scalar energy;
    scalar diameter;
};

// u(r) = -e / (1 + (r/d)^p)
struct softwell_potential<P=8> {
    scalar energy;
    scalar decay_distance;
};
```

### Threewise potentials

```c++
// u(r,s) = e dot(r,s) / |r||s|
struct bending_potential {
    scalar energy;
};
```


## Simulations

### Newtonian dynamics

```c++
struct newtonian_dynamics_config {
    scalar   timestep;
    step     steps;
    function callback;
};

void simulate_newtonian_dynamics(system, config);
```


### Brownian dynamics

```c++
struct brownian_dynamics_config {
    scalar   temperature;
    scalar   timestep;
    scalar   spacestep;
    step     steps;
    uint64_t seed;
    function callback;
};

void simulate_brownian_dynamics(system, config);
```
