//! Biodiversity metrics for the cultural ecosystem.
//!
//! Just as ecologists measure ecosystem health through diversity indices,
//! we measure the health of a cultural landscape. Higher biodiversity means
//! a more resilient cultural ecosystem — more perspectives, more innovation,
//! more capacity to adapt.

use serde::{Deserialize, Serialize};

use crate::TraditionSpecies;

/// A biodiversity report for a cultural ecosystem.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct BiodiversityReport {
    /// Shannon diversity index: H' = -Σ(pᵢ * ln(pᵢ))
    pub shannon_index: f64,
    /// Species (tradition) richness: count of non-extinct traditions.
    pub richness: usize,
    /// Pielou's evenness: J = H' / ln(S), where S = richness.
    pub evenness: f64,
}

impl BiodiversityReport {
    /// Compute biodiversity from a list of tradition populations.
    pub fn from_populations(populations: &[f64]) -> Self {
        let total: f64 = populations.iter().sum();
        if total <= 0.0 {
            return Self {
                shannon_index: 0.0,
                richness: 0,
                evenness: 0.0,
            };
        }

        let living: Vec<f64> = populations.iter().filter(|&&p| p > 0.0).cloned().collect();
        let richness = living.len();

        if richness == 0 {
            return Self {
                shannon_index: 0.0,
                richness: 0,
                evenness: 0.0,
            };
        }

        let shannon_index: f64 = living
            .iter()
            .map(|&p| {
                let pi = p / total;
                if pi > 0.0 { -pi * pi.ln() } else { 0.0 }
            })
            .sum();

        let evenness = if richness > 1 {
            shannon_index / (richness as f64).ln()
        } else {
            0.0
        };

        Self {
            shannon_index,
            richness,
            evenness,
        }
    }

    /// Compute biodiversity from a list of TraditionSpecies.
    pub fn from_species(species: &[TraditionSpecies], extinction_threshold: f64) -> Self {
        let pops: Vec<f64> = species
            .iter()
            .map(|s| if s.population > extinction_threshold { s.population } else { 0.0 })
            .collect();
        Self::from_populations(&pops)
    }

    /// Classify ecosystem health based on Shannon index.
    pub fn health(&self) -> EcosystemHealth {
        if self.richness == 0 {
            EcosystemHealth::Collapsed
        } else if self.shannon_index < 0.5 {
            EcosystemHealth::Depauperate
        } else if self.shannon_index < 1.0 {
            EcosystemHealth::Moderate
        } else if self.shannon_index < 2.0 {
            EcosystemHealth::Diverse
        } else {
            EcosystemHealth::Hyperdiverse
        }
    }

    /// Simpson's diversity index: D = 1 - Σ(pᵢ²)
    pub fn simpson_index(&self, populations: &[f64]) -> f64 {
        let total: f64 = populations.iter().sum();
        if total <= 0.0 {
            return 0.0;
        }
        1.0 - populations
            .iter()
            .map(|&p| {
                let pi = p / total;
                pi * pi
            })
            .sum::<f64>()
    }

    /// Berger-Parker dominance: proportion of the most abundant species.
    pub fn berger_parker(&self, populations: &[f64]) -> f64 {
        let total: f64 = populations.iter().sum();
        if total <= 0.0 {
            return 0.0;
        }
        let max_pop = populations.iter().cloned().fold(0.0_f64, f64::max);
        max_pop / total
    }
}

/// Ecosystem health classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EcosystemHealth {
    /// No living traditions.
    Collapsed,
    /// Very low diversity — monoculture or near-monoculture.
    Depauperate,
    /// Some diversity but dominated by few.
    Moderate,
    /// Good diversity, multiple coexisting traditions.
    Diverse,
    /// High diversity, many evenly distributed traditions.
    Hyperdiverse,
}

impl std::fmt::Display for EcosystemHealth {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EcosystemHealth::Collapsed => write!(f, "Collapsed"),
            EcosystemHealth::Depauperate => write!(f, "Depauperate"),
            EcosystemHealth::Moderate => write!(f, "Moderate"),
            EcosystemHealth::Diverse => write!(f, "Diverse"),
            EcosystemHealth::Hyperdiverse => write!(f, "Hyperdiverse"),
        }
    }
}

/// Compare biodiversity between two states (before/after perturbation).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct BiodiversityChange {
    pub initial: BiodiversityReport,
    pub final_: BiodiversityReport,
    pub shannon_delta: f64,
    pub richness_delta: i64,
    pub evenness_delta: f64,
    pub direction: DiversityTrend,
}

/// Direction of biodiversity change.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DiversityTrend {
    Increasing,
    Stable,
    Decreasing,
}

impl BiodiversityChange {
    pub fn compare(initial: BiodiversityReport, final_: BiodiversityReport) -> Self {
        let shannon_delta = final_.shannon_index - initial.shannon_index;
        let richness_delta = final_.richness as i64 - initial.richness as i64;
        let evenness_delta = final_.evenness - initial.evenness;

        let direction = if shannon_delta.abs() < 0.01 {
            DiversityTrend::Stable
        } else if shannon_delta > 0.0 {
            DiversityTrend::Increasing
        } else {
            DiversityTrend::Decreasing
        };

        Self {
            initial,
            final_,
            shannon_delta,
            richness_delta,
            evenness_delta,
            direction,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shannon_index_uniform() {
        let pops = vec![50.0, 50.0, 50.0, 50.0];
        let report = BiodiversityReport::from_populations(&pops);
        // Shannon for uniform: ln(4) ≈ 1.386
        assert!((report.shannon_index - 4.0_f64.ln()).abs() < 0.01);
        assert_eq!(report.richness, 4);
        assert!((report.evenness - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_shannon_index_monoculture() {
        let pops = vec![100.0, 0.0, 0.0];
        let report = BiodiversityReport::from_populations(&pops);
        assert!((report.shannon_index).abs() < 0.01);
        assert_eq!(report.richness, 1);
    }

    #[test]
    fn test_evenness_max_one() {
        let pops = vec![25.0, 25.0, 25.0, 25.0];
        let report = BiodiversityReport::from_populations(&pops);
        assert!((report.evenness - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_empty_ecosystem() {
        let pops: Vec<f64> = vec![];
        let report = BiodiversityReport::from_populations(&pops);
        assert_eq!(report.shannon_index, 0.0);
        assert_eq!(report.richness, 0);
    }

    #[test]
    fn test_all_zero_populations() {
        let pops = vec![0.0, 0.0, 0.0];
        let report = BiodiversityReport::from_populations(&pops);
        assert_eq!(report.richness, 0);
    }

    #[test]
    fn test_simpson_index() {
        let pops = vec![50.0, 50.0];
        let report = BiodiversityReport::from_populations(&pops);
        let simpson = report.simpson_index(&pops);
        // 1 - (0.25 + 0.25) = 0.5
        assert!((simpson - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_berger_parker() {
        let pops = vec![80.0, 20.0];
        let report = BiodiversityReport::from_populations(&pops);
        let bp = report.berger_parker(&pops);
        assert!((bp - 0.8).abs() < 1e-10);
    }

    #[test]
    fn test_ecosystem_health_levels() {
        assert_eq!(
            BiodiversityReport::from_populations(&[]).health(),
            EcosystemHealth::Collapsed
        );
        assert_eq!(
            BiodiversityReport::from_populations(&[99.0, 1.0]).health(),
            EcosystemHealth::Depauperate
        );
        assert_eq!(
            BiodiversityReport::from_populations(&[10.0; 10]).health(),
            EcosystemHealth::Hyperdiverse
        );
    }

    #[test]
    fn test_biodiversity_change_increasing() {
        let initial = BiodiversityReport::from_populations(&[99.0, 1.0]);
        let final_ = BiodiversityReport::from_populations(&[25.0, 25.0, 25.0, 25.0]);
        let change = BiodiversityChange::compare(initial, final_);
        assert_eq!(change.direction, DiversityTrend::Increasing);
        assert!(change.shannon_delta > 0.0);
    }

    #[test]
    fn test_biodiversity_change_decreasing() {
        let initial = BiodiversityReport::from_populations(&[25.0, 25.0, 25.0, 25.0]);
        let final_ = BiodiversityReport::from_populations(&[99.0, 1.0]);
        let change = BiodiversityChange::compare(initial, final_);
        assert_eq!(change.direction, DiversityTrend::Decreasing);
    }

    #[test]
    fn test_from_species() {
        let species = vec![
            TraditionSpecies::new("A", 50.0, 0.3, 100.0, 0.2),
            TraditionSpecies::new("B", 30.0, 0.3, 80.0, 0.8),
            TraditionSpecies::new("C", 0.1, 0.3, 50.0, 0.5), // effectively extinct
        ];
        let report = BiodiversityReport::from_species(&species, 1.0);
        assert_eq!(report.richness, 2);
    }

    #[test]
    fn test_health_display() {
        assert_eq!(EcosystemHealth::Hyperdiverse.to_string(), "Hyperdiverse");
    }
}
