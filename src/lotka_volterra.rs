use crate::error::EcologyError;
use crate::tradition::Tradition;

/// Solver method for the Lotka-Volterra ODE system.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SolverMethod {
    Euler,
    RK4,
}

use serde::{Deserialize, Serialize};

/// Configuration for a Lotka-Volterra competition simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LotkaVolterraConfig {
    pub traditions: Vec<Tradition>,
    pub competition_matrix: Vec<Vec<f64>>,
    pub dt: f64,
    pub method: SolverMethod,
}

impl LotkaVolterraConfig {
    /// Create a new config with validation.
    pub fn new(
        traditions: Vec<Tradition>,
        competition_matrix: Vec<Vec<f64>>,
        dt: f64,
        method: SolverMethod,
    ) -> Result<Self, EcologyError> {
        if traditions.is_empty() {
            return Err(EcologyError::EmptyTraditions);
        }
        if dt <= 0.0 {
            return Err(EcologyError::InvalidTimeStep { value: dt });
        }
        let n = traditions.len();
        if competition_matrix.len() != n
            || competition_matrix.iter().any(|row| row.len() != n)
        {
            return Err(EcologyError::MatrixDimensionMismatch {
                expected: n,
                got: competition_matrix.len(),
            });
        }
        // Validate dial dimensions are consistent
        let dim = traditions[0].dial_position.len();
        for t in &traditions {
            if t.dial_position.len() != dim {
                return Err(EcologyError::DialDimensionMismatch {
                    name: t.name.clone(),
                    expected: dim,
                    got: t.dial_position.len(),
                });
            }
            t.validate()?;
        }
        Ok(Self {
            traditions,
            competition_matrix,
            dt,
            method,
        })
    }

    /// Number of traditions (species).
    pub fn n(&self) -> usize {
        self.traditions.len()
    }

    /// Extract current population vector.
    pub fn populations(&self) -> Vec<f64> {
        self.traditions.iter().map(|t| t.population).collect()
    }
}

/// Result of a Lotka-Volterra simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationResult {
    /// Time series of population vectors: trajectory[step][tradition_idx]
    pub trajectory: Vec<Vec<f64>>,
    /// Time points
    pub times: Vec<f64>,
    /// Names of traditions in order
    pub tradition_names: Vec<String>,
}

/// Compute derivatives for the competitive Lotka-Volterra system:
/// dT_i/dt = r_i * T_i * (1 - Σ_j α_ij * T_j / K_i)
fn derivatives(populations: &[f64], config: &LotkaVolterraConfig) -> Vec<f64> {
    let n = populations.len();
    let mut derivs = vec![0.0; n];
    for i in 0..n {
        let r_i = config.traditions[i].growth_rate;
        let k_i = config.traditions[i].carrying_capacity;
        let competition_sum: f64 = (0..n)
            .map(|j| config.competition_matrix[i][j] * populations[j])
            .sum();
        derivs[i] = r_i * populations[i] * (1.0 - competition_sum / k_i);
    }
    derivs
}

/// Euler method single step.
fn euler_step(populations: &[f64], config: &LotkaVolterraConfig) -> Vec<f64> {
    let derivs = derivatives(populations, config);
    populations
        .iter()
        .zip(derivs.iter())
        .map(|(&p, &d)| (p + config.dt * d).max(0.0))
        .collect()
}

/// RK4 method single step.
fn rk4_step(populations: &[f64], config: &LotkaVolterraConfig) -> Vec<f64> {
    let dt = config.dt;
    let k1 = derivatives(populations, config);
    let p2: Vec<f64> = populations.iter().zip(k1.iter()).map(|(&p, &k)| p + 0.5 * dt * k).collect();
    let k2 = derivatives(&p2, config);
    let p3: Vec<f64> = populations.iter().zip(k2.iter()).map(|(&p, &k)| p + 0.5 * dt * k).collect();
    let k3 = derivatives(&p3, config);
    let p4: Vec<f64> = populations.iter().zip(k3.iter()).map(|(&p, &k)| p + dt * k).collect();
    let k4 = derivatives(&p4, config);

    populations
        .iter()
        .enumerate()
        .map(|(i, &p)| {
            (p + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])).max(0.0)
        })
        .collect()
}

/// Solve the Lotka-Volterra system for a given number of steps.
pub fn solve(config: &LotkaVolterraConfig, steps: usize) -> Result<SimulationResult, EcologyError> {
    let mut trajectory = Vec::with_capacity(steps + 1);
    let mut times = Vec::with_capacity(steps + 1);
    let tradition_names: Vec<String> = config.traditions.iter().map(|t| t.name.clone()).collect();

    let mut populations = config.populations();
    let mut t = 0.0;

    trajectory.push(populations.clone());
    times.push(t);

    for step in 1..=steps {
        populations = match config.method {
            SolverMethod::Euler => euler_step(&populations, config),
            SolverMethod::RK4 => rk4_step(&populations, config),
        };

        // Check for divergence
        if populations.iter().any(|&p| p > 1e10 || p.is_nan()) {
            return Err(EcologyError::Diverged { step });
        }

        t += config.dt;
        trajectory.push(populations.clone());
        times.push(t);
    }

    Ok(SimulationResult {
        trajectory,
        times,
        tradition_names,
    })
}

/// Solve until populations stabilize (max change < tolerance) or max_steps reached.
pub fn solve_to_equilibrium(
    config: &LotkaVolterraConfig,
    tolerance: f64,
    max_steps: usize,
) -> Result<SimulationResult, EcologyError> {
    let mut trajectory = vec![];
    let mut times = vec![];
    let tradition_names: Vec<String> = config.traditions.iter().map(|t| t.name.clone()).collect();

    let mut populations = config.populations();
    let mut t = 0.0;

    trajectory.push(populations.clone());
    times.push(t);

    for step in 1..=max_steps {
        let prev = populations.clone();
        populations = match config.method {
            SolverMethod::Euler => euler_step(&populations, config),
            SolverMethod::RK4 => rk4_step(&populations, config),
        };

        if populations.iter().any(|&p| p > 1e10 || p.is_nan()) {
            return Err(EcologyError::Diverged { step });
        }

        t += config.dt;
        trajectory.push(populations.clone());
        times.push(t);

        // Check convergence
        let max_change = populations
            .iter()
            .zip(prev.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        if max_change < tolerance {
            break;
        }
    }

    Ok(SimulationResult {
        trajectory,
        times,
        tradition_names,
    })
}
