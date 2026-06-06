//! Lotka-Volterra competition dynamics for traditions on the cultural dial.
//!
//! The classic two-species (or N-species) competition model:
//!   dT₁/dt = r₁T₁(1 - T₁/K₁ - α₁₂T₂/K₁)
//!   dT₂/dt = r₂T₂(1 - T₂/K₂ - α₂₁T₁/K₂)
//!
//! Traditions grow logistically but are suppressed by competitors proportional
//! to the competition coefficients (α) derived from niche overlap.

use serde::{Deserialize, Serialize};

use crate::TraditionSpecies;

/// Lotka-Volterra competition system for N traditions.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct LotkaVolterra {
    /// The competing tradition-species.
    pub species: Vec<TraditionSpecies>,
    /// N×N competition matrix. α[i][j] = effect of species j on species i.
    pub competition_matrix: Vec<Vec<f64>>,
}

impl LotkaVolterra {
    /// Create a new Lotka-Volterra system.
    pub fn new(
        species: Vec<TraditionSpecies>,
        competition_matrix: Vec<Vec<f64>>,
    ) -> Self {
        Self {
            species,
            competition_matrix,
        }
    }

    /// Create a system with only intraspecific competition (no intercompetition).
    /// All species grow independently to their carrying capacity.
    pub fn independent(species: Vec<TraditionSpecies>) -> Self {
        let n = species.len();
        let matrix = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| if i == j { 1.0 } else { 0.0 })
                    .collect()
            })
            .collect();
        Self::new(species, matrix)
    }

    /// Compute growth rates (dT/dt) for all species at current populations.
    pub fn growth_rates(&self) -> Vec<f64> {
        let n = self.species.len();
        let populations: Vec<f64> = self.species.iter().map(|s| s.population).collect();

        (0..n)
            .map(|i| {
                let si = &self.species[i];
                let mut competition_sum = 0.0;
                for j in 0..n {
                    competition_sum += self.competition_matrix[i][j] * populations[j];
                }
                si.growth_rate * si.population * (1.0 - competition_sum / si.carrying_capacity)
            })
            .collect()
    }

    /// Advance the system by one time step using Euler integration.
    pub fn step(&mut self, dt: f64) {
        let rates = self.growth_rates();
        for (i, s) in self.species.iter_mut().enumerate() {
            s.population = (s.population + rates[i] * dt).max(0.0);
        }
    }

    /// Simulate for `steps` iterations with time step `dt`.
    pub fn simulate(&mut self, steps: usize, dt: f64) -> Vec<Vec<f64>> {
        let mut history = Vec::with_capacity(steps + 1);
        history.push(self.species.iter().map(|s| s.population).collect());
        for _ in 0..steps {
            self.step(dt);
            history.push(self.species.iter().map(|s| s.population).collect());
        }
        history
    }

    /// Check if the system has reached equilibrium (all |dT/dt| < threshold).
    pub fn is_equilibrium(&self, threshold: f64) -> bool {
        self.growth_rates().iter().all(|&r| r.abs() < threshold)
    }

    /// Get total population across all species.
    pub fn total_population(&self) -> f64 {
        self.species.iter().map(|s| s.population).sum()
    }

    /// Find which species has the largest population.
    pub fn dominant_species(&self) -> Option<&TraditionSpecies> {
        self.species.iter().max_by(|a, b| {
            a.population
                .partial_cmp(&b.population)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Find species that have gone extinct (population < threshold).
    pub fn extinct_species(&self, threshold: f64) -> Vec<&TraditionSpecies> {
        self.species
            .iter()
            .filter(|s| s.population < threshold)
            .collect()
    }

    /// Compute the zero-growth isocline for species i.
    /// Returns (T_self, T_other) values where dT_i/dt = 0.
    /// For 2-species: T₁ = K₁ - α₁₂·T₂, T₂ = K₂ - α₂₁·T₁
    pub fn isocline_intercept(&self, species_idx: usize) -> (f64, f64) {
        let si = &self.species[species_idx];
        let alpha = &self.competition_matrix[species_idx];
        let mut interspecific_sum = 0.0;
        for (j, &a) in alpha.iter().enumerate() {
            if j != species_idx {
                interspecific_sum += a;
            }
        }
        // When all others at 0: T_self = K
        // When self at 0: sum of others = K / avg_alpha (simplified)
        (
            si.carrying_capacity,
            si.carrying_capacity / interspecific_sum.max(0.001),
        )
    }
}

/// Two-species Lotka-Volterra convenience function.
/// Returns growth rates (dT1/dt, dT2/dt).
pub fn two_species_rates(
    t1: &TraditionSpecies,
    t2: &TraditionSpecies,
    alpha12: f64,
    alpha21: f64,
) -> (f64, f64) {
    let dt1 = t1.growth_rate
        * t1.population
        * (1.0 - t1.population / t1.carrying_capacity - alpha12 * t2.population / t1.carrying_capacity);
    let dt2 = t2.growth_rate
        * t2.population
        * (1.0 - t2.population / t2.carrying_capacity - alpha21 * t1.population / t2.carrying_capacity);
    (dt1, dt2)
}

/// Simulate a simple two-species competition over time.
pub fn two_species_simulate(
    t1: TraditionSpecies,
    t2: TraditionSpecies,
    alpha12: f64,
    alpha21: f64,
    steps: usize,
    dt: f64,
) -> (Vec<f64>, Vec<f64>) {
    let mut t1 = t1;
    let mut t2 = t2;
    let mut pop1 = Vec::with_capacity(steps);
    let mut pop2 = Vec::with_capacity(steps);

    for _ in 0..steps {
        pop1.push(t1.population);
        pop2.push(t2.population);
        let (r1, r2) = two_species_rates(&t1, &t2, alpha12, alpha21);
        t1.population = (t1.population + r1 * dt).max(0.0);
        t2.population = (t2.population + r2 * dt).max(0.0);
    }

    (pop1, pop2)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_species(name: &str, pop: f64, r: f64, k: f64, pos: f64) -> TraditionSpecies {
        TraditionSpecies::new(name, pop, r, k, pos)
    }

    #[test]
    fn test_independent_growth_to_carrying_capacity() {
        let s = make_species("A", 10.0, 0.5, 100.0, 0.5);
        let mut lv = LotkaVolterra::independent(vec![s]);
        lv.simulate(1000, 0.1);
        let final_pop = lv.species[0].population;
        assert!((final_pop - 100.0).abs() < 1.0, "Expected ~100, got {final_pop}");
    }

    #[test]
    fn test_two_species_compete() {
        let t1 = make_species("A", 50.0, 0.3, 100.0, 0.3);
        let t2 = make_species("B", 50.0, 0.3, 100.0, 0.7);
        let matrix = vec![vec![1.0, 0.8], vec![0.8, 1.0]];
        let mut lv = LotkaVolterra::new(vec![t1, t2], matrix);
        lv.simulate(2000, 0.1);
        // Both should be below carrying capacity due to competition
        assert!(lv.species[0].population < 100.0);
        assert!(lv.species[1].population < 100.0);
    }

    #[test]
    fn test_growth_rates_logistic() {
        let s = make_species("A", 50.0, 0.5, 100.0, 0.5);
        let lv = LotkaVolterra::independent(vec![s]);
        let rates = lv.growth_rates();
        // dT/dt = r*T*(1 - T/K) = 0.5*50*(1 - 50/100) = 12.5
        assert!((rates[0] - 12.5).abs() < 1e-10);
    }

    #[test]
    fn test_at_carrying_capacity_rate_is_zero() {
        let s = make_species("A", 100.0, 0.5, 100.0, 0.5);
        let lv = LotkaVolterra::independent(vec![s]);
        let rates = lv.growth_rates();
        assert!(rates[0].abs() < 1e-10);
    }

    #[test]
    fn test_extinct_species_detection() {
        let t1 = make_species("A", 100.0, 0.3, 200.0, 0.3);
        let t2 = make_species("B", 0.1, 0.3, 200.0, 0.7);
        let lv = LotkaVolterra::independent(vec![t1, t2]);
        let extinct = lv.extinct_species(1.0);
        assert_eq!(extinct.len(), 1);
        assert_eq!(extinct[0].name, "B");
    }

    #[test]
    fn test_dominant_species() {
        let t1 = make_species("A", 50.0, 0.3, 100.0, 0.3);
        let t2 = make_species("B", 200.0, 0.3, 300.0, 0.7);
        let lv = LotkaVolterra::independent(vec![t1, t2]);
        assert_eq!(lv.dominant_species().unwrap().name, "B");
    }

    #[test]
    fn test_two_species_convenience_rates() {
        let t1 = make_species("A", 50.0, 0.5, 100.0, 0.3);
        let t2 = make_species("B", 30.0, 0.4, 80.0, 0.7);
        let (r1, _r2) = two_species_rates(&t1, &t2, 0.5, 0.5);
        // Should be negative competition effect
        let r1_no_comp = 0.5 * 50.0 * (1.0 - 50.0 / 100.0);
        assert!(r1 < r1_no_comp);
    }

    #[test]
    fn test_population_never_negative() {
        let t1 = make_species("A", 1.0, -0.5, 100.0, 0.5);
        let mut lv = LotkaVolterra::independent(vec![t1]);
        lv.simulate(100, 0.1);
        assert!(lv.species[0].population >= 0.0);
    }

    #[test]
    fn test_total_population() {
        let t1 = make_species("A", 50.0, 0.3, 100.0, 0.3);
        let t2 = make_species("B", 30.0, 0.3, 100.0, 0.7);
        let lv = LotkaVolterra::independent(vec![t1, t2]);
        assert!((lv.total_population() - 80.0).abs() < 1e-10);
    }

    #[test]
    fn test_equilibrium_detection() {
        let s = make_species("A", 100.0, 0.5, 100.0, 0.5);
        let lv = LotkaVolterra::independent(vec![s]);
        assert!(lv.is_equilibrium(0.01));
    }
}
