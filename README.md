# dial-ecology

**Ecological succession on the cultural dial.**

> *What if cultural traditions are species? What if the spectrum of human belief is a landscape — and every tradition is an organism competing for niche space, subject to the same merciless mathematics that governs forests, coral reefs, and tide pools?*

`dial-ecology` models cultural traditions as species in an ecological system. They grow, compete, coexist, and go extinct according to the laws of population ecology — specifically, the **Lotka-Volterra competition equations**. The "dial" is the spectrum of cultural expression: left to right, traditional to progressive, sacred to secular, old to new. Every tradition occupies a niche. Every niche has limited space.

## The Ecology of Culture

### 🌱 Traditions Are Species

A religion, a political ideology, a folk practice, a scientific paradigm — each is a "species" with:

- **Population** — number of adherents
- **Growth rate** — how fast it spreads (proselytizing religions grow fast; contemplative orders grow slow)
- **Carrying capacity** — the maximum population the niche can sustain
- **Dial position** — where on the cultural spectrum it lives

### ⚔️ Competition Is Real

Two traditions occupying the same niche **compete**. The Lotka-Volterra equations govern this:

```
dT₁/dt = r₁·T₁·(1 - T₁/K₁ - α₁₂·T₂/K₁)
dT₂/dt = r₂·T₂·(1 - T₂/K₂ - α₂₁·T₁/K₂)
```

Where:
- `r` = intrinsic growth rate
- `K` = carrying capacity
- `α` = competition coefficient (derived from niche overlap)

**High niche overlap → high α → fierce competition → one tradition may drive the other to extinction.**

### 🤝 Coexistence Requires Differentiation

Two traditions can coexist *if and only if*:

```
α₁₂ < K₁/K₂  AND  α₂₁ < K₂/K₁
```

This is the ecological translation of tolerance: *you can share the landscape, but only if you're different enough*. Identical niches mean war. Differentiated niches mean peace.

### 🌊 Succession Unfolds

Like bare rock colonized by lichens → mosses → grasses → shrubs → trees:

1. **Pioneer traditions** — fast-growing, low diversity, few adherents, high zeal
2. **Intermediate community** — more traditions emerge, competition intensifies
3. **Climax community** — stable, diverse, many coexisting traditions at equilibrium

The cultural landscape matures just like an ecosystem.

### 💀 Extinction Happens

A tradition goes extinct when its population falls below a critical threshold. The **Allee effect** makes this irreversible: small populations grow *slower*, not faster — they've lost the critical mass needed for cultural transmission. A language with 5 speakers is already dead. A religion with 50 practitioners is dying.

### 🌈 Biodiversity = Resilience

The **Shannon diversity index** measures ecosystem health:

```
H' = -Σ(pᵢ · ln(pᵢ))
```

- **High diversity** → many traditions, evenly distributed → resilient to disruption
- **Low diversity** → monoculture or near-monoculture → fragile, vulnerable to collapse
- **Monoculture** → one dominant tradition → maximally brittle

A healthy cultural ecosystem, like a healthy forest, is *diverse*.

## Modules

| Module | Description |
|--------|-------------|
| [`niche`](src/niche.rs) | Niche regions on the dial, overlap, competition coefficients, competitive exclusion |
| [`lotka_volterra`](src/lotka_volterra.rs) | Lotka-Volterra competition dynamics, simulation, equilibrium detection |
| [`succession`](src/succession.rs) | Ecological succession: pioneer → climax, stage classification |
| [`coexistence`](src/coexistence.rs) | Coexistence conditions, equilibrium analysis, pairwise compatibility |
| [`extinction`](src/extinction.rs) | Extinction dynamics, Allee effect, critical thresholds, risk assessment |
| [`biodiversity`](src/biodiversity.rs) | Shannon diversity, Simpson index, evenness, ecosystem health classification |

## Quick Start

```rust
use dial_ecology::{TraditionSpecies, lotka_volterra::LotkaVolterra, biodiversity::BiodiversityReport};

// Define competing traditions
let buddhism = TraditionSpecies::new("Buddhism", 500.0, 0.05, 1000.0, 0.4);
let materialism = TraditionSpecies::new("Materialism", 300.0, 0.08, 800.0, 0.7);

// Set up competition (moderate overlap)
let alpha_bv = 0.4; // Buddhism ← Materialism competition
let alpha_mb = 0.3; // Materialism ← Buddhism competition
let matrix = vec![vec![1.0, alpha_bv], vec![alpha_mb, 1.0]];

let mut ecosystem = LotkaVolterra::new(vec![buddhism, materialism], matrix);

// Simulate 1000 time steps
let history = ecosystem.simulate(1000, 0.1);

// Measure ecosystem health
let report = BiodiversityReport::from_species(&ecosystem.species, 1.0);
println!("Shannon diversity: {:.3}", report.shannon_index);
println!("Evenness: {:.3}", report.evenness);
println!("Ecosystem health: {}", report.health());
```

## Core Types

```rust
// A tradition modeled as an ecological species
struct TraditionSpecies {
    name: String,
    population: f64,        // current adherents
    growth_rate: f64,       // intrinsic rate of increase
    carrying_capacity: f64, // maximum sustainable population
    dial_position: f64,     // position on cultural spectrum [0, 1]
}

// The competition system
struct LotkaVolterra {
    species: Vec<TraditionSpecies>,
    competition_matrix: Vec<Vec<f64>>,  // α[i][j] = effect of j on i
}

// Succession stage
struct SuccessionStage {
    stage: usize,
    dominant: Vec<String>,
    diversity: f64,
}

// Coexistence analysis result
struct CoexistenceResult {
    can_coexist: bool,
    stable_equilibrium: Option<Vec<f64>>,
    winner: Option<usize>,
}

// Biodiversity metrics
struct BiodiversityReport {
    shannon_index: f64,  // H' = -Σ(pᵢ·ln(pᵢ))
    richness: usize,     // number of living traditions
    evenness: f64,       // J = H'/ln(S)
}
```

## Examples

### Competitive Exclusion

Two traditions in the *same niche* — one will win:

```rust
use dial_ecology::niche::{NicheRegion, competitive_exclusion};

let catholicism = NicheRegion::new(0.35, 0.1);
let protestantism = NicheRegion::new(0.37, 0.1);

// 80%+ overlap → competitive exclusion
assert!(competitive_exclusion(&catholicism, &protestantism));
```

### Peaceful Coexistence

Two traditions in *different niches* — both thrive:

```rust
use dial_ecology::coexistence::analyze_two_species;
use dial_ecology::TraditionSpecies;

let sufism = TraditionSpecies::new("Sufism", 100.0, 0.03, 200.0, 0.3);
let secular_humanism = TraditionSpecies::new("Secular Humanism", 400.0, 0.06, 600.0, 0.8);

let result = analyze_two_species(&sufism, &secular_humanism, 0.1, 0.1);
assert!(result.can_coexist); // Low competition → stable equilibrium
```

### Allee Effect and Extinction

A dying tradition spirals toward extinction:

```rust
use dial_ecology::extinction::{simulate_extinction, ExtinctionParams, ExtinctionRisk};
use dial_ecology::TraditionSpecies;

let dying_language = TraditionSpecies::new("Aramaic", 3.0, 0.02, 50.0, 0.2);
let params = ExtinctionParams::default();

// Assess risk
let risk = ExtinctionRisk::assess(&dying_language, &params);
assert_eq!(risk, ExtinctionRisk::High);

// Simulate
let result = simulate_extinction(dying_language, &params, 5000, 0.1);
// Small population → Allee effect → extinction spiral
```

## Installation

```toml
[dependencies]
dial-ecology = "0.1"
```

## Why?

Because culture is an ecosystem, and the math of ecology applies.

The same equations that describe how red squirrels outcompete grey squirrels in conifer forests also describe how Reformed theology displaced Catholicism in Scotland. The same Allee effect that dooms small populations of black-footed ferrets also dooms minority languages with 50 speakers. The same succession dynamics that turn lava flows into old-growth forests also turn frontier movements into established religions.

**Biodiversity isn't just biology. It's the measure of a healthy culture.**

A monoculture — one dominant tradition, no competition, no diversity — is *brittle*. It looks stable until it isn't. A diverse ecosystem — many traditions, competing and coexisting, each filling a niche — is *resilient*. It adapts. It survives shocks.

This library gives you the mathematical tools to think about culture ecologically. To quantify diversity. To predict when coexistence is possible and when it isn't. To model how traditions rise, compete, and fall — just like species in any ecosystem.

---

*In the ecology of the mind, as in the ecology of the forest, diversity is strength.*

## License

MIT
