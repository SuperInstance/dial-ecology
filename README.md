# dial-ecology

[![crates.io](https://img.shields.io/crates/v/dial-ecology.svg)](https://crates.io/crates/dial-ecology)
[![docs.rs](https://docs.rs/dial-ecology/badge.svg)](https://docs.rs/dial-ecology)
[![license: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## The Idea

Musical traditions compete for listener attention exactly the way species compete for ecological niches. Jazz and blues overlap in harmonic space but differ in rhythmic complexity. Electronic and classical occupy different spectral niches. When a new genre emerges, it either finds an empty niche (coexists), outcompetes a similar genre (replaces it), or can't establish (goes extinct).

The **Lotka-Volterra equations** model this precisely:

```
dTᵢ/dt = rᵢ · Tᵢ · (1 - Σⱼ αᵢⱼ · Tⱼ / K)
```

where Tᵢ is tradition i's popularity, rᵢ is its growth rate, K is carrying capacity (total audience), and αᵢⱼ is the competition coefficient — how much tradition j hinders tradition i. The competition coefficients come from **cultural dial positions**: traditions close on the dial compete hard; distant traditions barely interact.

## How It Works

### Define traditions on the cultural dial

```rust
use dial_ecology::Tradition;

let jazz = Tradition::new("Jazz", 0.8, 50.0)
    .with_dial_position(&[3.8, 3.5, 2.8]);  // harmonic tension, rhythmic complexity, spectral density

let electronic = Tradition::new("Electronic", 1.2, 60.0)
    .with_dial_position(&[2.0, 4.0, 4.5]);

let classical = Tradition::new("Classical", 0.5, 40.0)
    .with_dial_position(&[3.2, 2.0, 2.5]);
```

### Compute niche overlap

```rust
use dial_ecology::niche::NicheOverlap;

let overlap = NicheOverlap::from_traditions(&[jazz.clone(), electronic.clone(), classical.clone()]);
// αᵢⱼ = exp(-||dialᵢ - dialⱼ||² / σ²)
// Jazz ↔ Classical: high overlap (both high harmonic tension)
// Jazz ↔ Electronic: low overlap (different rhythmic + spectral profiles)
```

### Solve the dynamics

```rust
use dial_ecology::lotka_volterra::{LotkaVolterra, SolverMethod};

let mut lv = LotkaVolterra::new(vec![jazz, electronic, classical])
    .with_niche_overlap(&overlap)
    .dt(0.01)
    .method(SolverMethod::RK4);

// Integrate forward — watch traditions evolve
for _ in 0..5000 {
    lv.step();
}

let populations = lv.populations();
println!("Jazz: {:.1}, Electronic: {:.1}, Classical: {:.1}",
    populations[0], populations[1], populations[2]);
```

### Find equilibria

```rust
use dial_ecology::equilibrium::EquilibriumFinder;

let eq = EquilibriumFinder::find(&lv);
for e in &eq {
    println!("Equilibrium: {:?}", e.populations);
    println!("  Stable: {} (eigenvalues: {:?})", e.is_stable, e.eigenvalues);
}
```

Stability is determined by the Jacobian at the fixed point: all eigenvalues negative = stable attractor, any positive = unstable.

### Measure biodiversity

```rust
use dial_ecology::biodiversity::BiodiversityIndex;

let bio = BiodiversityIndex::from_populations(&populations);
println!("Shannon H = {:.3}", bio.shannon_index);      // higher = more diverse
println!("Simpson D = {:.3}", bio.simpson_index);      // probability two random listeners prefer same genre
println!("Pielou evenness = {:.3}", bio.evenness);      // 1 = perfectly balanced
```

Low Shannon H = one genre dominates (monoculture). High H = healthy ecosystem with coexistence.

### Succession modeling

```rust
use dial_ecology::succession::SuccessionModel;

let mut model = SuccessionModel::new(lv);
model.add_invasion("Hip-Hop", 10.0, 1.5, 55.0, &[3.0, 4.5, 3.0]); // new genre arrives
let history = model.simulate(1000); // 1000 time steps
// Watch: does Hip-Hop establish? Does it displace something? Coexist?
```

## Module Map

| Module | What it does |
|---|---|
| `tradition` | `Tradition` — musical tradition with dial position, growth rate, carrying capacity |
| `lotka_volterra` | `LotkaVolterra` — competitive LV system with Euler/RK4 integration |
| `niche` | `NicheOverlap` — competition coefficients from dial distance (Gaussian kernel) |
| `equilibrium` | `EquilibriumFinder` — fixed points + Jacobian stability analysis |
| `succession` | `SuccessionModel` — ecosystem evolution with invasion/extinction events |
| `biodiversity` | `BiodiversityIndex` — Shannon, Simpson, Berger-Parker, Pielou evenness |
| `error` | `EcologyError` |

## When To Use This

- **Music recommendation**: understand genre ecosystem dynamics, predict what's rising/falling
- **Cultural analysis**: model how traditions influence and compete with each other
- **Generative music**: use competition dynamics as a generative engine for evolving compositions
- **Agent fleet modeling**: agents competing for resources (compute, attention, bandwidth) follow the same dynamics

## Links

- [Documentation](https://docs.rs/dial-ecology)
- [Repository](https://github.com/SuperInstance/dial-ecology)
- [crates.io](https://crates.io/crates/dial-ecology)
- See also: [lotka-beats](https://crates.io/crates/lotka-beats) for LV dynamics in generative music

## License

MIT
