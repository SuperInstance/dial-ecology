use serde::{Deserialize, Serialize};

/// A musical tradition positioned on a cultural dial.
///
/// Each tradition has a population (popularity share), intrinsic growth rate,
/// carrying capacity, and a position vector on the cultural dial that
/// determines niche overlap with other traditions.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Tradition {
    pub name: String,
    pub population: f64,
    pub growth_rate: f64,
    pub carrying_capacity: f64,
    pub dial_position: Vec<f64>,
}

impl Tradition {
    /// Create a new tradition.
    pub fn new(
        name: impl Into<String>,
        population: f64,
        growth_rate: f64,
        carrying_capacity: f64,
        dial_position: Vec<f64>,
    ) -> Self {
        Self {
            name: name.into(),
            population,
            growth_rate,
            carrying_capacity,
            dial_position,
        }
    }

    /// Euclidean distance to another tradition on the cultural dial.
    pub fn distance_to(&self, other: &Tradition) -> f64 {
        self.dial_position
            .iter()
            .zip(other.dial_position.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt()
    }

    /// Validate this tradition's parameters.
    pub fn validate(&self) -> Result<(), crate::error::EcologyError> {
        if self.population < 0.0 {
            return Err(crate::error::EcologyError::NegativePopulation {
                name: self.name.clone(),
                value: self.population,
            });
        }
        if self.carrying_capacity <= 0.0 {
            return Err(crate::error::EcologyError::InvalidCarryingCapacity {
                name: self.name.clone(),
                value: self.carrying_capacity,
            });
        }
        Ok(())
    }
}
