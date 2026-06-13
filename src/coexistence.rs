//! Coexistence conditions for competing traditions.
//!
//! Two traditions can coexist when:
//!   α₁₂ < K₁/K₂  AND  α₂₁ < K₂/K₁
//!
//! When these conditions hold, the system has a stable interior equilibrium.
//! When they don't, one tradition drives the other to extinction (competitive exclusion).

use serde::{Deserialize, Serialize};

use crate::TraditionSpecies;
use crate::lotka_volterra::LotkaVolterra;

/// Result of a coexistence analysis.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct CoexistenceResult {
    /// Whether the species can coexist.
    pub can_coexist: bool,
    /// The stable equilibrium point, if one exists (populations at equilibrium).
    pub stable_equilibrium: Option<Vec<f64>>,
    /// Which species wins if they cannot coexist (index of winning species).
    pub winner: Option<usize>,
    /// Coexistence conditions checked: (alpha12_lt_K1_over_K2, alpha21_lt_K2_over_K1).
    pub conditions: Vec<bool>,
}

/// Analyze coexistence for a two-species Lotka-Volterra system.
pub fn analyze_two_species(
    t1: &TraditionSpecies,
    t2: &TraditionSpecies,
    alpha12: f64,
    alpha21: f64,
) -> CoexistenceResult {
    let k1 = t1.carrying_capacity;
    let k2 = t2.carrying_capacity;

    let cond1 = alpha12 < k1 / k2;
    let cond2 = alpha21 < k2 / k1;

    let can_coexist = cond1 && cond2;

    let (stable_equilibrium, winner) = if can_coexist {
        // Interior equilibrium: T₁* = (K₁ - α₁₂K₂) / (1 - α₁₂α₂₁)
        //                           T₂* = (K₂ - α₂₁K₁) / (1 - α₁₂α₂₁)
        let denom = 1.0 - alpha12 * alpha21;
        if denom.abs() > 1e-10 {
            let t1_eq = (k1 - alpha12 * k2) / denom;
            let t2_eq = (k2 - alpha21 * k1) / denom;
            if t1_eq > 0.0 && t2_eq > 0.0 {
                (Some(vec![t1_eq, t2_eq]), None)
            } else {
                // Both conditions met but equilibrium invalid — fall through
                let w = if alpha12 * k2 > k1 { 1 } else { 0 };
                (None, Some(w))
            }
        } else {
            (None, None)
        }
    } else {
        // Determine winner: species with higher relative carrying capacity wins
        let w = if !cond1 && !cond2 {
            // Both exceed — whoever has larger K wins
            if k1 >= k2 { 0 } else { 1 }
        } else if !cond1 {
            // alpha12 too high: species 2 wins
            1
        } else {
            // alpha21 too high: species 1 wins
            0
        };
        (None, Some(w))
    };

    CoexistenceResult {
        can_coexist,
        stable_equilibrium,
        winner,
        conditions: vec![cond1, cond2],
    }
}

/// Analyze coexistence for an N-species Lotka-Volterra system.
/// Uses simulation to detect stable equilibrium.
pub fn analyze_n_species(
    system: &LotkaVolterra,
    steps: usize,
    dt: f64,
    eq_threshold: f64,
) -> CoexistenceResult {
    let mut sys = system.clone();
    sys.simulate(steps, dt);

    let surviving: Vec<usize> = sys
        .species
        .iter()
        .enumerate()
        .filter(|(_, s)| s.population > eq_threshold)
        .map(|(i, _)| i)
        .collect();

    let can_coexist = surviving.len() > 1;
    let stable_equilibrium = if sys.is_equilibrium(eq_threshold * 10.0) {
        Some(sys.species.iter().map(|s| s.population).collect())
    } else {
        None
    };

    let winner = if !can_coexist && !surviving.is_empty() {
        Some(surviving[0])
    } else {
        None
    };

    CoexistenceResult {
        can_coexist,
        stable_equilibrium,
        winner,
        conditions: vec![sys.is_equilibrium(eq_threshold)],
    }
}

/// Check the Lotka-Volterra coexistence conditions for a competition matrix.
/// For species pair (i, j): αᵢⱼ < Kᵢ/Kⱼ AND αⱼᵢ < Kⱼ/Kᵢ
pub fn pairwise_coexistence_matrix(
    species: &[TraditionSpecies],
    competition_matrix: &[Vec<f64>],
) -> Vec<Vec<Option<bool>>> {
    let n = species.len();
    (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    if i == j {
                        None
                    } else {
                        let cond1 = competition_matrix[i][j] < species[i].carrying_capacity / species[j].carrying_capacity;
                        let cond2 = competition_matrix[j][i] < species[j].carrying_capacity / species[i].carrying_capacity;
                        Some(cond1 && cond2)
                    }
                })
                .collect()
        })
        .collect()
}

/// Compute the minimum niche differentiation required for coexistence.
/// Based on the competition coefficients: Δ_min = 1 - α₁₂·α₂₁ (simplified).
pub fn minimum_differentiation(alpha12: f64, alpha21: f64) -> f64 {
    (1.0 - alpha12 * alpha21).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sp(name: &str, pop: f64, r: f64, k: f64) -> TraditionSpecies {
        TraditionSpecies::new(name, pop, r, k, 0.5)
    }

    #[test]
    fn test_coexistence_weak_competition() {
        let t1 = sp("A", 50.0, 0.3, 100.0);
        let t2 = sp("B", 50.0, 0.3, 100.0);
        let result = analyze_two_species(&t1, &t2, 0.5, 0.5);
        assert!(result.can_coexist);
        assert!(result.stable_equilibrium.is_some());
        assert!(result.winner.is_none());
    }

    #[test]
    fn test_no_coexistence_strong_competition() {
        let t1 = sp("A", 50.0, 0.3, 100.0);
        let t2 = sp("B", 50.0, 0.3, 100.0);
        let result = analyze_two_species(&t1, &t2, 2.0, 2.0);
        assert!(!result.can_coexist);
        assert!(result.winner.is_some());
    }

    #[test]
    fn test_asymmetric_competition_species1_wins() {
        let t1 = sp("A", 50.0, 0.3, 100.0);
        let t2 = sp("B", 50.0, 0.3, 100.0);
        // alpha21 = 2.0 > K2/K1 = 1.0, so species 1 wins
        let result = analyze_two_species(&t1, &t2, 0.5, 2.0);
        assert!(!result.can_coexist);
        assert_eq!(result.winner, Some(0));
    }

    #[test]
    fn test_equilibrium_values_correct() {
        let t1 = sp("A", 10.0, 0.3, 100.0);
        let t2 = sp("B", 10.0, 0.3, 100.0);
        // α12 = 0.5, α21 = 0.5, K1 = K2 = 100
        // T1* = (100 - 0.5*100) / (1 - 0.25) = 50/0.75 = 66.67
        // T2* = same
        let result = analyze_two_species(&t1, &t2, 0.5, 0.5);
        let eq = result.stable_equilibrium.unwrap();
        assert!((eq[0] - 66.667).abs() < 0.1, "Expected ~66.67, got {}", eq[0]);
        assert!((eq[1] - 66.667).abs() < 0.1);
    }

    #[test]
    fn test_conditions_vector() {
        let t1 = sp("A", 50.0, 0.3, 100.0);
        let t2 = sp("B", 50.0, 0.3, 100.0);
        let result = analyze_two_species(&t1, &t2, 0.5, 0.5);
        assert_eq!(result.conditions, vec![true, true]);
    }

    #[test]
    fn test_n_species_coexistence() {
        let species = vec![
            TraditionSpecies::new("A", 50.0, 0.3, 100.0, 0.2),
            TraditionSpecies::new("B", 50.0, 0.3, 100.0, 0.8),
        ];
        let matrix = vec![vec![1.0, 0.3], vec![0.3, 1.0]];
        let system = LotkaVolterra::new(species, matrix);
        let result = analyze_n_species(&system, 5000, 0.1, 0.1);
        assert!(result.can_coexist);
    }

    #[test]
    fn test_pairwise_coexistence_matrix() {
        let species = vec![
            sp("A", 50.0, 0.3, 100.0),
            sp("B", 50.0, 0.3, 100.0),
        ];
        let matrix = vec![vec![1.0, 0.5], vec![0.5, 1.0]];
        let pw = pairwise_coexistence_matrix(&species, &matrix);
        assert_eq!(pw[0][1], Some(true));
        assert_eq!(pw[1][0], Some(true));
    }

    #[test]
    fn test_minimum_differentiation() {
        let diff = minimum_differentiation(0.3, 0.4);
        assert!((diff - 0.88).abs() < 1e-10);
    }

    #[test]
    fn test_minimum_differentiation_high_overlap() {
        let diff = minimum_differentiation(0.9, 0.9);
        assert!((diff - 0.19).abs() < 1e-10);
    }
}
