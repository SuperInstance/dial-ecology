use serde::{Deserialize, Serialize};

/// Biodiversity report for a genre ecosystem.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiodiversityReport {
    pub shannon_index: f64,
    pub simpson_index: f64,
    pub species_count: usize,
    pub evenness: f64,
}

/// Compute biodiversity indices from population counts.
///
/// Treats populations as absolute abundances and converts to proportions.
pub fn compute(populations: &[f64]) -> BiodiversityReport {
    let total: f64 = populations.iter().sum();
    if total <= 0.0 {
        return BiodiversityReport {
            shannon_index: 0.0,
            simpson_index: 0.0,
            species_count: 0,
            evenness: 0.0,
        };
    }

    let proportions: Vec<f64> = populations.iter().map(|&p| p / total).collect();

    // Shannon index: H = -Σ p_i * ln(p_i)
    let shannon = -proportions
        .iter()
        .filter(|&&p| p > 0.0)
        .map(|&p| p * p.ln())
        .sum::<f64>();

    // Simpson index: D = 1 - Σ p_i²
    let simpson = 1.0 - proportions.iter().map(|&p| p * p).sum::<f64>();

    // Species count (non-zero populations)
    let species_count = populations.iter().filter(|&&p| p > 0.0).count();

    // Pielou's evenness: J = H / ln(S)
    let evenness = if species_count > 1 {
        shannon / (species_count as f64).ln()
    } else {
        0.0
    };

    BiodiversityReport {
        shannon_index: shannon,
        simpson_index: simpson,
        species_count,
        evenness,
    }
}

/// Compute biodiversity over a trajectory (time series).
pub fn compute_trajectory(trajectory: &[Vec<f64>]) -> Vec<BiodiversityReport> {
    trajectory.iter().map(|pops| compute(pops)).collect()
}

/// Richness: count of traditions above a threshold.
pub fn richness(populations: &[f64], threshold: f64) -> usize {
    populations.iter().filter(|&&p| p > threshold).count()
}

/// Berger-Parker dominance index: proportion of the most abundant tradition.
pub fn berger_parker(populations: &[f64]) -> f64 {
    let total: f64 = populations.iter().sum();
    if total <= 0.0 {
        return 0.0;
    }
    let max_pop = populations.iter().cloned().fold(0.0_f64, f64::max);
    max_pop / total
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_uniform_distribution() {
        let pops = vec![100.0, 100.0, 100.0, 100.0];
        let report = compute(&pops);
        assert!((report.shannon_index - (4.0_f64).ln()).abs() < 1e-10);
        assert!((report.simpson_index - 0.75).abs() < 1e-10);
        assert!((report.evenness - 1.0).abs() < 1e-10);
        assert_eq!(report.species_count, 4);
    }

    #[test]
    fn test_single_dominant() {
        let pops = vec![1000.0, 1.0, 1.0];
        let report = compute(&pops);
        assert!(report.shannon_index < 0.5);
        assert!(report.simpson_index < 0.01);
        assert!(report.evenness < 0.5);
    }

    #[test]
    fn test_all_extinct() {
        let pops = vec![0.0, 0.0, 0.0];
        let report = compute(&pops);
        assert_eq!(report.shannon_index, 0.0);
        assert_eq!(report.species_count, 0);
    }
}
