use serde::{Deserialize, Serialize};

/// Niche overlap matrix computed from dial positions.
///
/// Traditions closer on the cultural dial have higher competition coefficients.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NicheOverlap {
    pub traditions: Vec<String>,
    pub overlap_matrix: Vec<Vec<f64>>,
}

/// Compute niche overlap (competition coefficients) from dial positions.
///
/// Uses Gaussian kernel: α_ij = exp(-d²_ij / (2σ²))
/// where d_ij is the Euclidean distance between dial positions.
///
/// Self-overlap (α_ii) is always 1.0.
pub fn compute_from_dial(
    names: Vec<String>,
    dial_positions: &[Vec<f64>],
    sigma: f64,
) -> NicheOverlap {
    let n = names.len();
    let sigma2 = sigma * sigma;
    let mut matrix = vec![vec![0.0; n]; n];

    for i in 0..n {
        matrix[i][i] = 1.0;
        for j in (i + 1)..n {
            let dist_sq: f64 = dial_positions[i]
                .iter()
                .zip(dial_positions[j].iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum();
            let overlap = (-dist_sq / (2.0 * sigma2)).exp();
            matrix[i][j] = overlap;
            matrix[j][i] = overlap;
        }
    }

    NicheOverlap {
        traditions: names,
        overlap_matrix: matrix,
    }
}

/// Compute competition coefficients from traditions directly.
pub fn compute_from_traditions(
    traditions: &[crate::tradition::Tradition],
    sigma: f64,
) -> NicheOverlap {
    let names: Vec<String> = traditions.iter().map(|t| t.name.clone()).collect();
    let positions: Vec<Vec<f64>> = traditions.iter().map(|t| t.dial_position.clone()).collect();
    compute_from_dial(names, &positions, sigma)
}

impl NicheOverlap {
    /// Get the competition coefficient between two traditions by name.
    pub fn get(&self, name_a: &str, name_b: &str) -> Option<f64> {
        let i = self.traditions.iter().position(|n| n == name_a)?;
        let j = self.traditions.iter().position(|n| n == name_b)?;
        Some(self.overlap_matrix[i][j])
    }

    /// Average niche overlap across all pairs (excluding self).
    pub fn mean_overlap(&self) -> f64 {
        let n = self.traditions.len();
        if n < 2 {
            return 0.0;
        }
        let mut sum = 0.0;
        let mut count = 0;
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    sum += self.overlap_matrix[i][j];
                    count += 1;
                }
            }
        }
        sum / count as f64
    }
}
