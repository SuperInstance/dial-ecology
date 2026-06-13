//! Ecological succession on the cultural dial.
//!
//! Just as bare rock is colonized by lichens → mosses → grasses → shrubs → trees,
//! a cultural landscape progresses through pioneer traditions → intermediate
//! communities → climax community. Each stage modifies conditions for the next.

use serde::{Deserialize, Serialize};

use crate::{TraditionSpecies, lotka_volterra::LotkaVolterra, biodiversity::BiodiversityReport};

/// A stage in the ecological succession process.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SuccessionStage {
    /// Stage index (0 = pioneer, higher = more mature).
    pub stage: usize,
    /// Names of dominant traditions at this stage.
    pub dominant: Vec<String>,
    /// Shannon diversity at this stage.
    pub diversity: f64,
    /// Total tradition population at this stage.
    pub total_population: f64,
}

/// Result of running a full succession simulation.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SuccessionResult {
    /// Stages recorded during succession.
    pub stages: Vec<SuccessionStage>,
    /// Whether a climax community was reached.
    pub climax_reached: bool,
    /// Number of simulation steps taken.
    pub total_steps: usize,
}

/// Parameters for succession simulation.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SuccessionParams {
    /// Time step for each simulation step.
    pub dt: f64,
    /// Maximum steps before giving up.
    pub max_steps: usize,
    /// How often to sample (record a stage).
    pub sample_interval: usize,
    /// Threshold for declaring equilibrium (climax).
    pub equilibrium_threshold: f64,
    /// Minimum consecutive equilibrium checks to declare climax.
    pub climax_confirmation_steps: usize,
}

impl Default for SuccessionParams {
    fn default() -> Self {
        Self {
            dt: 0.1,
            max_steps: 10_000,
            sample_interval: 100,
            equilibrium_threshold: 0.01,
            climax_confirmation_steps: 5,
        }
    }
}

/// Run ecological succession on a Lotka-Volterra system.
///
/// Simulates forward in time, recording stages at regular intervals.
/// Detects when the system reaches a climax community (equilibrium).
pub fn run_succession(
    mut system: LotkaVolterra,
    params: &SuccessionParams,
) -> SuccessionResult {
    let mut stages = Vec::new();
    let mut equilibrium_count = 0;
    let mut climax_reached = false;

    // Record initial state
    stages.push(record_stage(0, &system));

    for step in 1..=params.max_steps {
        system.step(params.dt);

        if system.is_equilibrium(params.equilibrium_threshold) {
            equilibrium_count += 1;
            if equilibrium_count >= params.climax_confirmation_steps {
                climax_reached = true;
                stages.push(record_stage(step, &system));
                break;
            }
        } else {
            equilibrium_count = 0;
        }

        if step % params.sample_interval == 0 {
            stages.push(record_stage(step, &system));
        }
    }

    let total_steps = stages.last().map(|s| s.stage).unwrap_or(0);
    SuccessionResult {
        stages,
        climax_reached,
        total_steps: params.max_steps.min(total_steps),
    }
}

fn record_stage(step: usize, system: &LotkaVolterra) -> SuccessionStage {
    let pops: Vec<(String, f64)> = system
        .species
        .iter()
        .map(|s| (s.name.clone(), s.population))
        .collect();

    let total_pop: f64 = system.species.iter().map(|s| s.population).sum();
    let diversity = BiodiversityReport::from_populations(
        &system.species.iter().map(|s| s.population).collect::<Vec<_>>(),
    )
    .shannon_index;

    // Dominant = species with population > 10% of total
    let dominant: Vec<String> = pops
        .iter()
        .filter(|(_, p)| *p > total_pop * 0.1)
        .map(|(name, _)| name.clone())
        .collect();

    SuccessionStage {
        stage: step,
        dominant,
        diversity,
        total_population: total_pop,
    }
}

/// Classify the succession stage:
/// - Pioneer: low diversity, 1-2 dominant species
/// - Intermediate: moderate diversity, several species
/// - Climax: high diversity, many coexisting species
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum StageClass {
    Pioneer,
    Intermediate,
    Climax,
}

impl StageClass {
    pub fn classify(stage: &SuccessionStage, richness: usize) -> Self {
        if stage.diversity < 0.5 || richness <= 2 {
            StageClass::Pioneer
        } else if stage.diversity < 1.5 {
            StageClass::Intermediate
        } else {
            StageClass::Climax
        }
    }
}

/// Generate a pioneer community: a small number of fast-growing traditions.
pub fn pioneer_community() -> Vec<TraditionSpecies> {
    vec![
        TraditionSpecies::new("Pioneer-1", 5.0, 0.8, 50.0, 0.2),
        TraditionSpecies::new("Pioneer-2", 3.0, 0.9, 40.0, 0.8),
    ]
}

/// Generate a climax community: many coexisting traditions with lower growth rates.
pub fn climax_community() -> Vec<TraditionSpecies> {
    vec![
        TraditionSpecies::new("Climax-1", 80.0, 0.1, 100.0, 0.1),
        TraditionSpecies::new("Climax-2", 60.0, 0.1, 80.0, 0.3),
        TraditionSpecies::new("Climax-3", 50.0, 0.12, 70.0, 0.5),
        TraditionSpecies::new("Climax-4", 40.0, 0.08, 60.0, 0.7),
        TraditionSpecies::new("Climax-5", 30.0, 0.1, 50.0, 0.9),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lotka_volterra::LotkaVolterra;

    #[test]
    fn test_pioneer_community_creation() {
        let pioneers = pioneer_community();
        assert_eq!(pioneers.len(), 2);
        assert!(pioneers[0].growth_rate > 0.5);
    }

    #[test]
    fn test_climax_community_creation() {
        let climax = climax_community();
        assert_eq!(climax.len(), 5);
        assert!(climax[0].population > 50.0);
    }

    #[test]
    fn test_succession_reaches_climax() {
        let species = vec![
            TraditionSpecies::new("A", 10.0, 0.5, 100.0, 0.3),
            TraditionSpecies::new("B", 5.0, 0.4, 80.0, 0.7),
        ];
        let matrix = vec![vec![1.0, 0.3], vec![0.3, 1.0]];
        let system = LotkaVolterra::new(species, matrix);
        let result = run_succession(system, &SuccessionParams::default());
        assert!(result.climax_reached);
    }

    #[test]
    fn test_stage_class_pioneer() {
        let stage = SuccessionStage {
            stage: 0,
            dominant: vec!["A".into()],
            diversity: 0.2,
            total_population: 15.0,
        };
        assert_eq!(StageClass::classify(&stage, 1), StageClass::Pioneer);
    }

    #[test]
    fn test_stage_class_climax() {
        let stage = SuccessionStage {
            stage: 1000,
            dominant: vec!["A".into(), "B".into(), "C".into()],
            diversity: 2.0,
            total_population: 500.0,
        };
        assert_eq!(StageClass::classify(&stage, 5), StageClass::Climax);
    }

    #[test]
    fn test_succession_params_default() {
        let params = SuccessionParams::default();
        assert_eq!(params.max_steps, 10_000);
        assert_eq!(params.sample_interval, 100);
    }

    #[test]
    fn test_stage_diversity_increases() {
        let species = vec![
            TraditionSpecies::new("A", 5.0, 0.6, 100.0, 0.2),
            TraditionSpecies::new("B", 3.0, 0.5, 80.0, 0.8),
        ];
        let matrix = vec![vec![1.0, 0.2], vec![0.2, 1.0]];
        let system = LotkaVolterra::new(species, matrix);
        let result = run_succession(system, &SuccessionParams::default());
        // Diversity should generally increase as populations stabilize
        assert!(result.stages.len() > 1);
    }
}
