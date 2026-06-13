//! Niche space on the cultural dial.
//!
//! Each tradition occupies a region on the dial [0, 1]. Niche overlap between
//! two traditions determines competition strength — the more overlap, the
//! fiercer the competition. Implements the competitive exclusion principle:
//! no two traditions can occupy the exact same niche indefinitely.

use serde::{Deserialize, Serialize};

use crate::{TraditionSpecies, DIAL_MIN, DIAL_MAX};

/// A niche region on the dial, defined by center position and width.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct NicheRegion {
    /// Center of the niche on the dial [0, 1].
    pub center: f64,
    /// Half-width of the niche; total span is [center - width, center + width].
    pub width: f64,
}

impl NicheRegion {
    /// Create a new niche region, clamping to dial bounds.
    pub fn new(center: f64, width: f64) -> Self {
        let center = center.clamp(DIAL_MIN, DIAL_MAX);
        let width = width.max(0.0);
        Self { center, width }
    }

    /// Left edge of the niche on the dial.
    pub fn left(&self) -> f64 {
        (self.center - self.width).max(DIAL_MIN)
    }

    /// Right edge of the niche on the dial.
    pub fn right(&self) -> f64 {
        (self.center + self.width).min(DIAL_MAX)
    }

    /// Total niche breadth (span).
    pub fn breadth(&self) -> f64 {
        self.right() - self.left()
    }

    /// Check if a dial position falls within this niche.
    pub fn contains(&self, position: f64) -> bool {
        position >= self.left() && position <= self.right()
    }

    /// Overlap area between two niche regions.
    pub fn overlap(&self, other: &NicheRegion) -> f64 {
        let left = self.left().max(other.left());
        let right = self.right().min(other.right());
        (right - left).max(0.0)
    }

    /// Fraction of this niche that overlaps with another (competition strength).
    /// Returns 0.0 if no overlap, up to 1.0 if completely overlapped.
    pub fn overlap_fraction(&self, other: &NicheRegion) -> f64 {
        if self.breadth() == 0.0 {
            return 0.0;
        }
        self.overlap(other) / self.breadth()
    }

    /// Competition coefficient derived from niche overlap.
    /// Uses Pianka's index: α = (overlap) / sqrt(breadth_self * breadth_other).
    pub fn competition_coefficient(&self, other: &NicheRegion) -> f64 {
        let denom = (self.breadth() * other.breadth()).sqrt();
        if denom == 0.0 {
            return 0.0;
        }
        self.overlap(other) / denom
    }
}

/// Compute the niche region from a tradition's dial position with a given width.
pub fn tradition_niche(tradition: &TraditionSpecies, niche_width: f64) -> NicheRegion {
    NicheRegion::new(tradition.dial_position, niche_width)
}

/// Build a full competition matrix from a set of traditions and their niche widths.
/// Returns an N×N matrix where entry (i, j) is the competition coefficient αᵢⱼ.
pub fn competition_matrix(traditions: &[TraditionSpecies], niche_widths: &[f64]) -> Vec<Vec<f64>> {
    let n = traditions.len();
    let niches: Vec<NicheRegion> = traditions
        .iter()
        .zip(niche_widths.iter())
        .map(|(t, &w)| tradition_niche(t, w))
        .collect();

    (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    if i == j {
                        1.0 // Intraspecific competition
                    } else {
                        niches[i].competition_coefficient(&niches[j])
                    }
                })
                .collect()
        })
        .collect()
}

/// Competitive exclusion check: if two traditions have >90% niche overlap,
/// they cannot coexist indefinitely.
pub fn competitive_exclusion(niche_a: &NicheRegion, niche_b: &NicheRegion) -> bool {
    let overlap = niche_a.overlap_fraction(niche_b);
    let reverse_overlap = niche_b.overlap_fraction(niche_a);
    overlap > 0.9 || reverse_overlap > 0.9
}

/// Niche differentiation: minimum distance between niche centers for coexistence.
/// Returns the minimum separation needed given the niche widths.
pub fn minimum_separation(width_a: f64, width_b: f64) -> f64 {
    (width_a + width_b) * 0.5
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::TraditionSpecies;

    #[test]
    fn test_niche_region_new_clamps_center() {
        let n = NicheRegion::new(1.5, 0.2);
        assert_eq!(n.center, 1.0);
        let n2 = NicheRegion::new(-0.5, 0.2);
        assert_eq!(n2.center, 0.0);
    }

    #[test]
    fn test_niche_region_edges() {
        let n = NicheRegion::new(0.5, 0.2);
        assert!((n.left() - 0.3).abs() < 1e-10);
        assert!((n.right() - 0.7).abs() < 1e-10);
    }

    #[test]
    fn test_niche_breadth() {
        let n = NicheRegion::new(0.5, 0.3);
        assert!((n.breadth() - 0.6).abs() < 1e-10);
    }

    #[test]
    fn test_niche_contains() {
        let n = NicheRegion::new(0.5, 0.2);
        assert!(n.contains(0.4));
        assert!(n.contains(0.5));
        assert!(n.contains(0.7));
        assert!(!n.contains(0.71));
    }

    #[test]
    fn test_niche_overlap_partial() {
        let a = NicheRegion::new(0.3, 0.2); // [0.1, 0.5]
        let b = NicheRegion::new(0.5, 0.2); // [0.3, 0.7]
        assert!((a.overlap(&b) - 0.2).abs() < 1e-10);
    }

    #[test]
    fn test_niche_overlap_none() {
        let a = NicheRegion::new(0.2, 0.1); // [0.1, 0.3]
        let b = NicheRegion::new(0.8, 0.1); // [0.7, 0.9]
        assert_eq!(a.overlap(&b), 0.0);
    }

    #[test]
    fn test_overlap_fraction() {
        let a = NicheRegion::new(0.5, 0.2); // breadth 0.4
        let b = NicheRegion::new(0.5, 0.2); // breadth 0.4
        // Full overlap
        assert!((a.overlap_fraction(&b) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_competition_coefficient_symmetric() {
        let a = NicheRegion::new(0.4, 0.2); // [0.2, 0.6], breadth 0.4
        let b = NicheRegion::new(0.6, 0.2); // [0.4, 0.8], breadth 0.4
        let ab = a.competition_coefficient(&b);
        let ba = b.competition_coefficient(&a);
        assert!((ab - ba).abs() < 1e-10);
    }

    #[test]
    fn test_competitive_exclusion_identical() {
        let a = NicheRegion::new(0.5, 0.2);
        let b = NicheRegion::new(0.5, 0.2);
        assert!(competitive_exclusion(&a, &b));
    }

    #[test]
    fn test_competitive_exclusion_separated() {
        let a = NicheRegion::new(0.2, 0.05);
        let b = NicheRegion::new(0.8, 0.05);
        assert!(!competitive_exclusion(&a, &b));
    }

    #[test]
    fn test_competition_matrix_diagonal_is_one() {
        let t = TraditionSpecies::new("A", 100.0, 0.1, 1000.0, 0.5);
        let traditions = vec![t];
        let widths = vec![0.2];
        let m = competition_matrix(&traditions, &widths);
        assert_eq!(m[0][0], 1.0);
    }

    #[test]
    fn test_minimum_separation() {
        let sep = minimum_separation(0.2, 0.3);
        assert!((sep - 0.25).abs() < 1e-10);
    }

    #[test]
    fn test_niche_clamps_to_dial_bounds() {
        let n = NicheRegion::new(0.05, 0.2);
        assert!((n.left() - 0.0).abs() < 1e-10);
        let n2 = NicheRegion::new(0.95, 0.2);
        assert!((n2.right() - 1.0).abs() < 1e-10);
    }
}
