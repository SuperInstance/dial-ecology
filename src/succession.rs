use crate::lotka_volterra::{self, LotkaVolterraConfig, SolverMethod};
use crate::niche::{self, NicheOverlap};
use crate::tradition::Tradition;
use serde::{Deserialize, Serialize};

/// A snapshot of the genre ecosystem at a point in time.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EcosystemSnapshot {
    pub time: f64,
    pub populations: Vec<f64>,
    pub tradition_names: Vec<String>,
}

/// Result of a full succession simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SuccessionResult {
    pub snapshots: Vec<EcosystemSnapshot>,
    pub initial_populations: Vec<f64>,
    pub final_populations: Vec<f64>,
    pub dominant_tradition: String,
    pub extinction_events: Vec<ExtinctionEvent>,
    pub niche_overlap: NicheOverlap,
}

/// Record of a tradition going (near-)extinct.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExtinctionEvent {
    pub tradition: String,
    pub time: f64,
    pub final_population: f64,
}

/// Succession model: simulate genre ecosystem evolution.
pub struct SuccessionModel {
    pub traditions: Vec<Tradition>,
    pub sigma: f64,
    pub dt: f64,
    pub method: SolverMethod,
}

impl SuccessionModel {
    pub fn new(
        traditions: Vec<Tradition>,
        sigma: f64,
        dt: f64,
        method: SolverMethod,
    ) -> Self {
        Self { traditions, sigma, dt, method }
    }

    /// Run the full succession simulation for `steps` time steps.
    pub fn simulate(&self, steps: usize, extinction_threshold: f64) -> SuccessionResult {
        let niche = niche::compute_from_traditions(&self.traditions, self.sigma);
        let config = LotkaVolterraConfig {
            traditions: self.traditions.clone(),
            competition_matrix: niche.overlap_matrix.clone(),
            dt: self.dt,
            method: self.method,
        };

        let result = lotka_volterra::solve(&config, steps).expect("simulation should not fail");

        let initial_populations = result.trajectory.first().unwrap().clone();
        let final_populations = result.trajectory.last().unwrap().clone();

        // Find dominant tradition at end
        let dominant_idx = final_populations
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        let dominant_tradition = result.tradition_names[dominant_idx].clone();

        // Detect extinction events
        let mut extinctions = Vec::new();
        let n = result.tradition_names.len();
        let mut went_extinct = vec![false; n];

        for (step_idx, populations) in result.trajectory.iter().enumerate() {
            for i in 0..n {
                if !went_extinct[i] && populations[i] < extinction_threshold && initial_populations[i] > extinction_threshold {
                    went_extinct[i] = true;
                    extinctions.push(ExtinctionEvent {
                        tradition: result.tradition_names[i].clone(),
                        time: result.times[step_idx],
                        final_population: populations[i],
                    });
                }
            }
        }

        let snapshots: Vec<EcosystemSnapshot> = result
            .trajectory
            .iter()
            .zip(result.times.iter())
            .map(|(pops, &t)| EcosystemSnapshot {
                time: t,
                populations: pops.clone(),
                tradition_names: result.tradition_names.clone(),
            })
            .collect();

        SuccessionResult {
            snapshots,
            initial_populations,
            final_populations,
            dominant_tradition,
            extinction_events: extinctions,
            niche_overlap: niche,
        }
    }

    /// Add a new tradition (invasion) and simulate from current state.
    pub fn simulate_invasion(
        &self,
        invader: Tradition,
        current_populations: &[f64],
        steps: usize,
    ) -> SuccessionResult {
        let mut all_traditions = self.traditions.clone();
        // Update existing populations
        for (i, t) in all_traditions.iter_mut().enumerate() {
            if i < current_populations.len() {
                t.population = current_populations[i];
            }
        }
        all_traditions.push(invader);

        let model = SuccessionModel::new(
            all_traditions,
            self.sigma,
            self.dt,
            self.method,
        );
        model.simulate(steps, 0.01)
    }
}
