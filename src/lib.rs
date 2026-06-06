//! # dial-ecology
//!
//! Ecological succession on the cultural dial: traditions compete like species
//! for niche space, governed by Lotka-Volterra dynamics.
//!
//! This library models cultural traditions as species in an ecological system,
//! where they compete for "space" on the cultural dial — the spectrum of human
//! expression, belief, and practice. Just as biological species compete for
//! resources, traditions compete for adherents, attention, and legitimacy.

use serde::{Deserialize, Serialize};

pub mod biodiversity;
pub mod coexistence;
pub mod extinction;
pub mod lotka_volterra;
pub mod niche;
pub mod succession;

/// Minimum position on the dial.
pub const DIAL_MIN: f64 = 0.0;
/// Maximum position on the dial.
pub const DIAL_MAX: f64 = 1.0;

/// A tradition modeled as an ecological species.
///
/// Each tradition has a population (number of adherents), a growth rate,
/// a carrying capacity (maximum sustainable population), and a position
/// on the cultural dial [0, 1].
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TraditionSpecies {
    /// Name of the tradition.
    pub name: String,
    /// Current population (number of adherents).
    pub population: f64,
    /// Intrinsic growth rate (r in Lotka-Volterra).
    pub growth_rate: f64,
    /// Carrying capacity (K) — maximum sustainable population.
    pub carrying_capacity: f64,
    /// Position on the cultural dial [0, 1].
    pub dial_position: f64,
}

impl TraditionSpecies {
    /// Create a new tradition species.
    pub fn new(
        name: impl Into<String>,
        population: f64,
        growth_rate: f64,
        carrying_capacity: f64,
        dial_position: f64,
    ) -> Self {
        Self {
            name: name.into(),
            population: population.max(0.0),
            growth_rate,
            carrying_capacity: carrying_capacity.max(0.0),
            dial_position: dial_position.clamp(DIAL_MIN, DIAL_MAX),
        }
    }

    /// Proportion of carrying capacity currently used.
    pub fn capacity_usage(&self) -> f64 {
        if self.carrying_capacity == 0.0 {
            return 0.0;
        }
        self.population / self.carrying_capacity
    }

    /// Check if the tradition is effectively extinct (below threshold).
    pub fn is_extinct(&self, threshold: f64) -> bool {
        self.population < threshold
    }

    /// Population as a fraction of carrying capacity [0, 1+].
    pub fn relative_abundance(&self) -> f64 {
        if self.carrying_capacity == 0.0 {
            return 0.0;
        }
        self.population / self.carrying_capacity
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tradition_new() {
        let t = TraditionSpecies::new("Test", 100.0, 0.5, 200.0, 0.5);
        assert_eq!(t.name, "Test");
        assert_eq!(t.population, 100.0);
        assert_eq!(t.growth_rate, 0.5);
        assert_eq!(t.carrying_capacity, 200.0);
        assert_eq!(t.dial_position, 0.5);
    }

    #[test]
    fn test_tradition_clamps_dial_position() {
        let t = TraditionSpecies::new("T", 10.0, 0.1, 100.0, 1.5);
        assert_eq!(t.dial_position, 1.0);
        let t2 = TraditionSpecies::new("T", 10.0, 0.1, 100.0, -0.5);
        assert_eq!(t2.dial_position, 0.0);
    }

    #[test]
    fn test_tradition_negative_population_clamped() {
        let t = TraditionSpecies::new("T", -10.0, 0.1, 100.0, 0.5);
        assert_eq!(t.population, 0.0);
    }

    #[test]
    fn test_capacity_usage() {
        let t = TraditionSpecies::new("T", 50.0, 0.1, 100.0, 0.5);
        assert!((t.capacity_usage() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_is_extinct() {
        let t = TraditionSpecies::new("T", 0.5, 0.1, 100.0, 0.5);
        assert!(t.is_extinct(1.0));
        assert!(!t.is_extinct(0.1));
    }

    #[test]
    fn test_relative_abundance() {
        let t = TraditionSpecies::new("T", 75.0, 0.1, 100.0, 0.5);
        assert!((t.relative_abundance() - 0.75).abs() < 1e-10);
    }
}
