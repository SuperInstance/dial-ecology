//! # dial-ecology
//!
//! Lotka-Volterra dynamics for musical tradition competition on cultural dials.
//!
//! Musical traditions compete for listener attention like species compete for
//! resources. This crate models tradition population (popularity) dynamics using
//! competitive Lotka-Volterra equations:
//!
//! `dT_i/dt = r_i · T_i · (1 - Σ_j α_ij · T_j / K_i)`
//!
//! where:
//! - `T_i` is the population (popularity) of tradition *i*
//! - `r_i` is the intrinsic growth rate
//! - `K_i` is the carrying capacity
//! - `α_ij` is the competition coefficient (niche overlap)
//!
//! ## Modules
//!
//! - [`tradition`] — Musical traditions with dial positions
//! - [`lotka_volterra`] — ODE solver (Euler, RK4)
//! - [`niche`] — Niche overlap from cultural dial distances
//! - [`equilibrium`] — Fixed point and stability analysis
//! - [`succession`] — Ecosystem evolution simulation
//! - [`biodiversity`] — Shannon, Simpson diversity indices

pub mod biodiversity;
pub mod error;
pub mod equilibrium;
pub mod lotka_volterra;
pub mod niche;
pub mod succession;
pub mod tradition;

pub use biodiversity::{BiodiversityReport, compute as compute_biodiversity};
pub use error::EcologyError;
pub use equilibrium::Equilibrium;
pub use lotka_volterra::{LotkaVolterraConfig, SimulationResult, SolverMethod, solve, solve_to_equilibrium};
pub use niche::NicheOverlap;
pub use succession::{SuccessionModel, SuccessionResult};
pub use tradition::Tradition;
