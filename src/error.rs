use std::fmt;

/// Errors that can occur in ecological simulations.
#[derive(Debug)]
pub enum EcologyError {
    EmptyTraditions,
    MatrixDimensionMismatch { expected: usize, got: usize },
    DialDimensionMismatch { name: String, expected: usize, got: usize },
    InvalidGrowthRate { name: String, value: f64 },
    NegativePopulation { name: String, value: f64 },
    InvalidCarryingCapacity { name: String, value: f64 },
    InvalidTimeStep { value: f64 },
    Diverged { step: usize },
    NoEquilibriumFound { iterations: usize },
}

impl fmt::Display for EcologyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EcologyError::EmptyTraditions => write!(f, "tradition list is empty"),
            EcologyError::MatrixDimensionMismatch { expected, got } => {
                write!(f, "competition matrix dimension mismatch: expected {expected}x{expected}, got {got}x{got}")
            }
            EcologyError::DialDimensionMismatch { name, expected, got } => {
                write!(f, "dial position dimension mismatch for tradition '{name}': expected {expected}, got {got}")
            }
            EcologyError::InvalidGrowthRate { name, value } => {
                write!(f, "invalid growth rate for tradition '{name}': {value}")
            }
            EcologyError::NegativePopulation { name, value } => {
                write!(f, "negative population for tradition '{name}': {value}")
            }
            EcologyError::InvalidCarryingCapacity { name, value } => {
                write!(f, "non-positive carrying capacity for tradition '{name}': {value}")
            }
            EcologyError::InvalidTimeStep { value } => {
                write!(f, "time step dt must be positive, got {value}")
            }
            EcologyError::Diverged { step } => {
                write!(f, "solver diverged: population exceeded bounds at step {step}")
            }
            EcologyError::NoEquilibriumFound { iterations } => {
                write!(f, "no equilibrium found after {iterations} iterations")
            }
        }
    }
}

impl std::error::Error for EcologyError {}
