//! Extinction dynamics for traditions on the cultural dial.
//!
//! When does a tradition die out? Several factors drive cultural extinction:
//! - Population falling below a critical threshold
//! - Allee effect: small populations have negative growth (can't sustain themselves)
//! - Competitive exclusion by a dominant tradition
//! - Environmental change shifting the cultural landscape

use serde::{Deserialize, Serialize};

use crate::TraditionSpecies;

/// Parameters controlling extinction dynamics.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ExtinctionParams {
    /// Population below which a tradition is considered extinct.
    pub extinction_threshold: f64,
    /// Allee effect threshold: populations below this have reduced growth.
    pub allee_threshold: f64,
    /// Critical minimum population for recovery (Allee effect strength).
    pub allee_strength: f64,
}

impl Default for ExtinctionParams {
    fn default() -> Self {
        Self {
            extinction_threshold: 1.0,
            allee_threshold: 10.0,
            allee_strength: 0.5,
        }
    }
}

/// Result of an extinction analysis.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ExtinctionResult {
    /// Whether the tradition went extinct during simulation.
    pub went_extinct: bool,
    /// Time step at which extinction occurred (if applicable).
    pub extinction_step: Option<usize>,
    /// Minimum population reached during simulation.
    pub min_population: f64,
    /// Whether the tradition experienced an Allee effect.
    pub allee_effect_detected: bool,
    /// Final population.
    pub final_population: f64,
}

/// Growth rate with Allee effect.
/// Below the Allee threshold, growth rate is reduced (or becomes negative).
/// dT/dt = r * T * (1 - T/K) * (T / (T + A)) where A is Allee threshold.
pub fn growth_with_allee(
    tradition: &TraditionSpecies,
    params: &ExtinctionParams,
) -> f64 {
    let t = tradition.population;
    let r = tradition.growth_rate;
    let k = tradition.carrying_capacity;

    let logistic = r * t * (1.0 - t / k);
    let allee_factor = t / (t + params.allee_threshold);

    logistic * allee_factor
}

/// Growth rate with strong Allee effect (negative growth below threshold).
/// dT/dt = r * T * (1 - T/K) * (T/A - 1) where A is critical threshold.
pub fn growth_with_strong_allee(
    tradition: &TraditionSpecies,
    critical_threshold: f64,
) -> f64 {
    let t = tradition.population;
    let r = tradition.growth_rate;
    let k = tradition.carrying_capacity;

    let logistic = r * t * (1.0 - t / k);
    let allee_factor = (t / critical_threshold) - 1.0;

    logistic * allee_factor
}

/// Simulate a single tradition with Allee effect, returning extinction trajectory.
pub fn simulate_extinction(
    mut tradition: TraditionSpecies,
    params: &ExtinctionParams,
    steps: usize,
    dt: f64,
) -> ExtinctionResult {
    let mut min_population = tradition.population;
    let mut extinction_step = None;
    let mut allee_effect_detected = false;

    for step in 0..steps {
        if tradition.population < params.extinction_threshold {
            extinction_step = Some(step);
            break;
        }

        if tradition.population < params.allee_threshold {
            allee_effect_detected = true;
        }

        let rate = growth_with_allee(&tradition, params);
        tradition.population = (tradition.population + rate * dt).max(0.0);
        min_population = min_population.min(tradition.population);
    }

    let went_extinct = tradition.population < params.extinction_threshold
        || extinction_step.is_some();

    ExtinctionResult {
        went_extinct,
        extinction_step,
        min_population,
        allee_effect_detected,
        final_population: tradition.population,
    }
}

/// Compute the critical threshold for a tradition: the population below which
/// recovery is impossible given competition from other traditions.
pub fn critical_threshold(
    tradition: &TraditionSpecies,
    competitor_population: f64,
    alpha: f64,
) -> f64 {
    // At critical threshold, growth rate = 0 despite competition
    // r*T*(1 - T/K - α*C/K) = 0, T ≠ 0
    // 1 - T/K - α*C/K = 0
    // T = K - α*C
    let threshold = tradition.carrying_capacity - alpha * competitor_population;
    threshold.max(0.0)
}

/// Time to extinction estimate for a declining tradition.
/// Uses exponential decay approximation: T(t) ≈ T₀ * e^(r_eff * t)
/// where r_eff is the effective (negative) growth rate.
pub fn time_to_extinction(
    current_population: f64,
    effective_growth_rate: f64,
    threshold: f64,
) -> Option<f64> {
    if effective_growth_rate >= 0.0 || current_population <= threshold {
        return None; // Won't go extinct or already extinct
    }
    // T₀ * e^(r*t) = threshold => t = ln(threshold/T₀) / r
    if current_population <= 0.0 || threshold <= 0.0 {
        return None;
    }
    let t = (threshold / current_population).ln() / effective_growth_rate;
    Some(t.max(0.0))
}

/// Assess extinction risk for a tradition based on current state.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ExtinctionRisk {
    /// Population healthy, no immediate risk.
    Low,
    /// Population declining or small.
    Moderate,
    /// Population very small, Allee effect likely.
    High,
    /// Population at or below threshold.
    Critical,
}

impl ExtinctionRisk {
    pub fn assess(tradition: &TraditionSpecies, params: &ExtinctionParams) -> Self {
        let pop = tradition.population;
        let threshold = params.extinction_threshold;
        let allee = params.allee_threshold;

        if pop <= threshold {
            ExtinctionRisk::Critical
        } else if pop < allee * 0.5 {
            ExtinctionRisk::High
        } else if pop < allee || tradition.growth_rate < 0.0 {
            ExtinctionRisk::Moderate
        } else {
            ExtinctionRisk::Low
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sp(name: &str, pop: f64, r: f64, k: f64) -> TraditionSpecies {
        TraditionSpecies::new(name, pop, r, k, 0.5)
    }

    #[test]
    fn test_growth_with_allee_normal() {
        let t = sp("A", 50.0, 0.5, 100.0);
        let params = ExtinctionParams::default();
        let rate = growth_with_allee(&t, &params);
        // Should be positive for a healthy population
        assert!(rate > 0.0);
    }

    #[test]
    fn test_growth_with_allee_small_population() {
        let t = sp("A", 2.0, 0.5, 100.0);
        let params = ExtinctionParams::default();
        let rate = growth_with_allee(&t, &params);
        // Allee factor = 2/(2+10) = 0.167, should be reduced
        let normal_rate = 0.5 * 2.0 * (1.0 - 2.0 / 100.0);
        assert!(rate < normal_rate);
    }

    #[test]
    fn test_strong_allee_below_threshold() {
        let t = sp("A", 5.0, 0.5, 100.0);
        let rate = growth_with_strong_allee(&t, 10.0);
        // Below critical: allee_factor = 5/10 - 1 = -0.5, negative growth
        assert!(rate < 0.0);
    }

    #[test]
    fn test_strong_allee_above_threshold() {
        let t = sp("A", 50.0, 0.5, 100.0);
        let rate = growth_with_strong_allee(&t, 10.0);
        assert!(rate > 0.0);
    }

    #[test]
    fn test_simulate_healthy_no_extinction() {
        let t = sp("A", 50.0, 0.3, 100.0);
        let result = simulate_extinction(t, &ExtinctionParams::default(), 1000, 0.1);
        assert!(!result.went_extinct);
    }

    #[test]
    fn test_simulate_tiny_population_extinction() {
        let t = sp("A", 0.5, 0.1, 100.0);
        let mut params = ExtinctionParams::default();
        params.allee_threshold = 5.0;
        params.extinction_threshold = 0.1;
        // Use strong Allee effect to ensure extinction
        let result = simulate_extinction_strong(t, &params, 10.0, 10000, 0.1);
        assert!(result.went_extinct);
    }

    fn simulate_extinction_strong(
        mut tradition: TraditionSpecies,
        params: &ExtinctionParams,
        critical_threshold: f64,
        steps: usize,
        dt: f64,
    ) -> ExtinctionResult {
        let mut min_population = tradition.population;
        let mut extinction_step = None;
        let mut allee_effect_detected = false;

        for step in 0..steps {
            if tradition.population < params.extinction_threshold {
                extinction_step = Some(step);
                break;
            }

            if tradition.population < params.allee_threshold {
                allee_effect_detected = true;
            }

            let rate = growth_with_strong_allee(&tradition, critical_threshold);
            tradition.population = (tradition.population + rate * dt).max(0.0);
            min_population = min_population.min(tradition.population);
        }

        let went_extinct = tradition.population < params.extinction_threshold
            || extinction_step.is_some();

        ExtinctionResult {
            went_extinct,
            extinction_step,
            min_population,
            allee_effect_detected,
            final_population: tradition.population,
        }
    }

    #[test]
    fn test_critical_threshold_no_competition() {
        let t = sp("A", 50.0, 0.3, 100.0);
        let thresh = critical_threshold(&t, 0.0, 0.0);
        assert_eq!(thresh, 100.0); // K - 0 = K
    }

    #[test]
    fn test_critical_threshold_with_competition() {
        let t = sp("A", 50.0, 0.3, 100.0);
        let thresh = critical_threshold(&t, 50.0, 0.5);
        // K - α*C = 100 - 25 = 75
        assert_eq!(thresh, 75.0);
    }

    #[test]
    fn test_time_to_extinction_declining() {
        let t = time_to_extinction(100.0, -0.1, 1.0);
        assert!(t.is_some());
        let t = t.unwrap();
        // ln(1/100) / -0.1 = ln(0.01) / -0.1 ≈ 46.05
        assert!(t > 45.0 && t < 47.0);
    }

    #[test]
    fn test_time_to_extinction_growing_is_none() {
        let t = time_to_extinction(100.0, 0.1, 1.0);
        assert!(t.is_none());
    }

    #[test]
    fn test_extinction_risk_levels() {
        let params = ExtinctionParams::default();
        assert_eq!(
            ExtinctionRisk::assess(&sp("A", 0.5, 0.3, 100.0), &params),
            ExtinctionRisk::Critical
        );
        assert_eq!(
            ExtinctionRisk::assess(&sp("A", 3.0, 0.3, 100.0), &params),
            ExtinctionRisk::High
        );
        assert_eq!(
            ExtinctionRisk::assess(&sp("A", 8.0, 0.3, 100.0), &params),
            ExtinctionRisk::Moderate
        );
        assert_eq!(
            ExtinctionRisk::assess(&sp("A", 50.0, 0.3, 100.0), &params),
            ExtinctionRisk::Low
        );
    }
}
