use dial_ecology::*;

fn make_tradition(name: &str, pop: f64, growth: f64, k: f64, pos: Vec<f64>) -> Tradition {
    Tradition::new(name, pop, growth, k, pos)
}

fn two_far_traditions() -> Vec<Tradition> {
    vec![
        make_tradition("Jazz", 50.0, 0.5, 100.0, vec![0.0, 0.0]),
        make_tradition("Metal", 50.0, 0.5, 100.0, vec![10.0, 10.0]),
    ]
}

fn two_close_traditions() -> Vec<Tradition> {
    vec![
        make_tradition("Jazz", 50.0, 0.5, 100.0, vec![0.0, 0.0]),
        make_tradition("Blues", 50.0, 0.5, 100.0, vec![0.1, 0.1]),
    ]
}

fn three_traditions() -> Vec<Tradition> {
    vec![
        make_tradition("Jazz", 50.0, 0.5, 100.0, vec![0.0, 0.0]),
        make_tradition("Blues", 40.0, 0.4, 100.0, vec![1.0, 0.0]),
        make_tradition("Rock", 30.0, 0.6, 100.0, vec![5.0, 5.0]),
    ]
}

fn make_config(traditions: Vec<Tradition>, method: SolverMethod) -> LotkaVolterraConfig {
    let niche = niche::compute_from_traditions(&traditions, 1.0);
    LotkaVolterraConfig::new(traditions, niche.overlap_matrix, 0.01, method).unwrap()
}

// ===== Tradition tests =====

#[test]
fn tradition_new() {
    let t = make_tradition("Test", 10.0, 0.5, 100.0, vec![1.0, 2.0]);
    assert_eq!(t.name, "Test");
    assert_eq!(t.population, 10.0);
    assert_eq!(t.growth_rate, 0.5);
    assert_eq!(t.carrying_capacity, 100.0);
    assert_eq!(t.dial_position, vec![1.0, 2.0]);
}

#[test]
fn tradition_distance_same() {
    let t = make_tradition("A", 10.0, 0.5, 100.0, vec![1.0, 2.0]);
    assert_eq!(t.distance_to(&t), 0.0);
}

#[test]
fn tradition_distance_different() {
    let a = make_tradition("A", 10.0, 0.5, 100.0, vec![0.0, 0.0]);
    let b = make_tradition("B", 10.0, 0.5, 100.0, vec![3.0, 4.0]);
    assert!((a.distance_to(&b) - 5.0).abs() < 1e-10);
}

#[test]
fn tradition_validate_negative_population() {
    let t = make_tradition("Bad", -1.0, 0.5, 100.0, vec![0.0]);
    assert!(t.validate().is_err());
}

#[test]
fn tradition_validate_zero_capacity() {
    let t = make_tradition("Bad", 10.0, 0.5, 0.0, vec![0.0]);
    assert!(t.validate().is_err());
}

#[test]
fn tradition_validate_ok() {
    let t = make_tradition("Good", 10.0, 0.5, 100.0, vec![0.0]);
    assert!(t.validate().is_ok());
}

#[test]
fn tradition_serialization() {
    let t = make_tradition("Jazz", 50.0, 0.5, 100.0, vec![1.0, 2.0]);
    let json = serde_json::to_string(&t).unwrap();
    let t2: Tradition = serde_json::from_str(&json).unwrap();
    assert_eq!(t.name, t2.name);
    assert_eq!(t.population, t2.population);
    assert_eq!(t.dial_position, t2.dial_position);
}

// ===== Config validation tests =====

#[test]
fn config_empty_traditions() {
    let result = LotkaVolterraConfig::new(vec![], vec![], 0.01, SolverMethod::Euler);
    assert!(result.is_err());
}

#[test]
fn config_bad_dt() {
    let traditions = two_far_traditions();
    let result = LotkaVolterraConfig::new(traditions, vec![vec![1.0, 0.0], vec![0.0, 1.0]], -0.1, SolverMethod::Euler);
    assert!(result.is_err());
}

#[test]
fn config_matrix_mismatch() {
    let traditions = two_far_traditions();
    let result = LotkaVolterraConfig::new(traditions, vec![vec![1.0]], 0.01, SolverMethod::Euler);
    assert!(result.is_err());
}

#[test]
fn config_dial_dimension_mismatch() {
    let traditions = vec![
        make_tradition("A", 10.0, 0.5, 100.0, vec![0.0, 0.0]),
        make_tradition("B", 10.0, 0.5, 100.0, vec![0.0]),  // different dimension
    ];
    let result = LotkaVolterraConfig::new(
        traditions,
        vec![vec![1.0, 0.0], vec![0.0, 1.0]],
        0.01,
        SolverMethod::Euler,
    );
    assert!(result.is_err());
}

#[test]
fn config_valid() {
    let config = make_config(two_far_traditions(), SolverMethod::Euler);
    assert_eq!(config.n(), 2);
}

// ===== Niche overlap tests =====

#[test]
fn niche_self_overlap_is_one() {
    let niche = niche::compute_from_traditions(&two_far_traditions(), 1.0);
    for i in 0..niche.overlap_matrix.len() {
        assert!((niche.overlap_matrix[i][i] - 1.0).abs() < 1e-10);
    }
}

#[test]
fn niche_far_traditions_low_overlap() {
    let niche = niche::compute_from_traditions(&two_far_traditions(), 1.0);
    let overlap = niche.overlap_matrix[0][1];
    assert!(overlap < 0.001, "far traditions should have near-zero overlap, got {overlap}");
}

#[test]
fn niche_close_traditions_high_overlap() {
    let niche = niche::compute_from_traditions(&two_close_traditions(), 1.0);
    let overlap = niche.overlap_matrix[0][1];
    assert!(overlap > 0.9, "close traditions should have high overlap, got {overlap}");
}

#[test]
fn niche_symmetric() {
    let niche = niche::compute_from_traditions(&three_traditions(), 1.0);
    let n = niche.overlap_matrix.len();
    for i in 0..n {
        for j in 0..n {
            assert!((niche.overlap_matrix[i][j] - niche.overlap_matrix[j][i]).abs() < 1e-10);
        }
    }
}

#[test]
fn niche_get_by_name() {
    let niche = niche::compute_from_traditions(&two_far_traditions(), 1.0);
    let val = niche.get("Jazz", "Metal");
    assert!(val.is_some());
    assert!(val.unwrap() < 0.001);
    assert!(niche.get("Jazz", "Nonexistent").is_none());
}

#[test]
fn niche_mean_overlap() {
    let niche = niche::compute_from_traditions(&two_close_traditions(), 1.0);
    let mean = niche.mean_overlap();
    assert!(mean > 0.9);
}

// ===== Solver tests =====

#[test]
fn solve_euler_basic() {
    let config = make_config(two_far_traditions(), SolverMethod::Euler);
    let result = solve(&config, 100).unwrap();
    assert_eq!(result.trajectory.len(), 101);
    assert_eq!(result.times.len(), 101);
    assert_eq!(result.tradition_names, vec!["Jazz", "Metal"]);
}

#[test]
fn solve_rk4_basic() {
    let config = make_config(two_far_traditions(), SolverMethod::RK4);
    let result = solve(&config, 100).unwrap();
    assert_eq!(result.trajectory.len(), 101);
}

#[test]
fn solve_populations_non_negative() {
    let config = make_config(three_traditions(), SolverMethod::RK4);
    let result = solve(&config, 1000).unwrap();
    for pops in &result.trajectory {
        for &p in pops {
            assert!(p >= 0.0, "population went negative: {p}");
        }
    }
}

#[test]
fn solve_no_competition_reaches_capacity() {
    // With identity competition matrix and low initial pop, should grow toward K
    let traditions = vec![
        make_tradition("Jazz", 10.0, 1.0, 100.0, vec![0.0]),
        make_tradition("Metal", 10.0, 1.0, 100.0, vec![100.0]),
    ];
    let config = LotkaVolterraConfig::new(
        traditions,
        vec![vec![1.0, 0.0], vec![0.0, 1.0]], // no interspecific competition
        0.01,
        SolverMethod::RK4,
    ).unwrap();
    let result = solve(&config, 5000).unwrap();
    let final_pops = result.trajectory.last().unwrap();
    for &p in final_pops {
        assert!((p - 100.0).abs() < 1.0, "should reach carrying capacity, got {p}");
    }
}

#[test]
fn solve_to_equilibrium_converges() {
    let config = make_config(two_far_traditions(), SolverMethod::RK4);
    let result = solve_to_equilibrium(&config, 1e-6, 100000).unwrap();
    // Should converge before max_steps
    assert!(result.trajectory.len() < 100000);
}

// ===== Competitive exclusion tests =====

#[test]
fn competitive_exclusion_close_traditions() {
    // Two very close traditions with different fitness: stronger should dominate
    let traditions = vec![
        make_tradition("Jazz", 50.0, 1.0, 100.0, vec![0.0, 0.0]),
        make_tradition("Blues", 50.0, 0.5, 100.0, vec![0.1, 0.1]),
    ];
    // High competition (close together) with asymmetric growth
    let config = LotkaVolterraConfig::new(
        traditions,
        vec![vec![1.0, 0.99], vec![0.99, 1.0]],
        0.01,
        SolverMethod::RK4,
    ).unwrap();
    let result = solve(&config, 100000).unwrap();
    let final_pops = result.trajectory.last().unwrap();
    // Stronger grower (Jazz, r=1.0) should dominate weaker (Blues, r=0.5)
    assert!(final_pops[0] > final_pops[1],
        "expected Jazz to dominate Blues, got {:?}", final_pops);
}

#[test]
fn coexistence_far_traditions() {
    // Far apart traditions should coexist
    let config = make_config(two_far_traditions(), SolverMethod::RK4);
    let result = solve(&config, 50000).unwrap();
    let final_pops = result.trajectory.last().unwrap();
    for &p in final_pops {
        assert!(p > 10.0, "tradition went extinct: {p}");
    }
}

// ===== Equilibrium tests =====

#[test]
fn equilibrium_stable_coexistence() {
    // No interspecific competition → stable coexistence at carrying capacities
    let traditions = vec![
        make_tradition("Jazz", 50.0, 0.5, 100.0, vec![0.0]),
        make_tradition("Metal", 50.0, 0.5, 100.0, vec![100.0]),
    ];
    let config = LotkaVolterraConfig::new(
        traditions,
        vec![vec![1.0, 0.0], vec![0.0, 1.0]],
        0.01,
        SolverMethod::RK4,
    ).unwrap();
    let eq = equilibrium::find_equilibrium(&config, 1e-6, 100000).unwrap();
    assert!(eq.is_stable);
    assert!((eq.populations[0] - 100.0).abs() < 1.0);
    assert!((eq.populations[1] - 100.0).abs() < 1.0);
}

#[test]
fn boundary_equilibria_found() {
    let config = make_config(two_far_traditions(), SolverMethod::Euler);
    let eqs = equilibrium::find_boundary_equilibria(&config);
    // Should have: trivial (all zero) + one per tradition
    assert_eq!(eqs.len(), 3);
    // Trivial equilibrium
    assert!(eqs[0].populations.iter().all(|&p| p == 0.0));
}

// ===== Biodiversity tests =====

#[test]
fn biodiversity_uniform_high() {
    let report = biodiversity::compute(&[25.0, 25.0, 25.0, 25.0]);
    assert!(report.shannon_index > 1.0);
    assert!(report.evenness > 0.99);
    assert_eq!(report.species_count, 4);
}

#[test]
fn biodiversity_dominant_low() {
    let report = biodiversity::compute(&[1000.0, 1.0, 1.0]);
    assert!(report.shannon_index < 0.5);
    assert!(report.simpson_index < 0.01);
}

#[test]
fn biodiversity_evenness_perfect() {
    let report = biodiversity::compute(&[10.0, 10.0, 10.0]);
    assert!((report.evenness - 1.0).abs() < 1e-10);
}

#[test]
fn biodiversity_trajectory() {
    let reports = biodiversity::compute_trajectory(&[
        vec![10.0, 10.0],
        vec![20.0, 20.0],
        vec![50.0, 50.0],
    ]);
    assert_eq!(reports.len(), 3);
    assert!((reports[0].shannon_index - reports[1].shannon_index).abs() < 1e-10);
}

#[test]
fn biodiversity_richness() {
    assert_eq!(biodiversity::richness(&[10.0, 0.0, 5.0], 0.1), 2);
    assert_eq!(biodiversity::richness(&[0.0, 0.0], 0.1), 0);
}

#[test]
fn biodiversity_berger_parker() {
    let bp = biodiversity::berger_parker(&[75.0, 25.0]);
    assert!((bp - 0.75).abs() < 1e-10);
}

// ===== Succession tests =====

#[test]
fn succession_basic() {
    let model = succession::SuccessionModel::new(three_traditions(), 1.0, 0.01, SolverMethod::RK4);
    let result = model.simulate(1000, 0.01);
    assert_eq!(result.snapshots.len(), 1001);
    assert!(!result.dominant_tradition.is_empty());
}

#[test]
fn succession_extinction_detection() {
    // Asymmetric competition: strong competitor suppresses weak one
    let traditions = vec![
        make_tradition("Dominant", 50.0, 1.0, 100.0, vec![0.0]),
        make_tradition("Weak", 50.0, 0.5, 100.0, vec![0.0]),
    ];
    // Dominant barely affected by Weak, but Weak strongly affected by Dominant
    let config = LotkaVolterraConfig::new(
        traditions,
        vec![vec![1.0, 0.2], vec![1.5, 1.0]],
        0.01,
        SolverMethod::RK4,
    ).unwrap();
    let result = solve(&config, 100000).unwrap();
    let final_pops = result.trajectory.last().unwrap();
    assert!(final_pops[0] > final_pops[1] * 3.0,
        "dominant should be much bigger: {:?}", final_pops);
}

#[test]
fn succession_invasion() {
    let model = succession::SuccessionModel::new(
        vec![make_tradition("Jazz", 80.0, 0.5, 100.0, vec![0.0])],
        1.0,
        0.01,
        SolverMethod::RK4,
    );
    let invader = make_tradition("HipHop", 5.0, 0.8, 100.0, vec![5.0, 5.0]);
    let result = model.simulate_invasion(invader, &[80.0], 5000);
    assert_eq!(result.snapshots[0].tradition_names.len(), 2);
    assert!(result.snapshots[0].populations.contains(&5.0));
}

// ===== Integration tests =====

#[test]
fn oscillations_possible() {
    // Asymmetric competition with overshoot: tradition A starts near capacity,
    // B starts low. B grows fast, overshoots, then settles.
    // With high interspecific competition on A, A drops below equilibrium first.
    let traditions = vec![
        make_tradition("A", 90.0, 1.0, 100.0, vec![0.0]),
        make_tradition("B", 10.0, 3.0, 100.0, vec![0.0]),
    ];
    let config = LotkaVolterraConfig::new(
        traditions,
        vec![vec![1.0, 0.5], vec![2.0, 1.0]], // asymmetric: B competes strongly with A
        0.001,
        SolverMethod::RK4,
    ).unwrap();
    let result = solve(&config, 200000).unwrap();

    // Tradition B starts at 10, should overshoot its equilibrium then come back
    let b_pops: Vec<f64> = result.trajectory.iter().map(|p| p[1]).collect();
    let final_b = b_pops.last().unwrap();
    let max_b = b_pops.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    // B should overshoot its final equilibrium value
    assert!(max_b > final_b + 1.0,
        "expected B to overshoot equilibrium: max={max_b}, final={final_b}");
}

#[test]
fn full_pipeline() {
    // End-to-end: create traditions → compute niche → solve → find equilibrium → biodiversity
    let traditions = three_traditions();
    let niche = niche::compute_from_traditions(&traditions, 2.0);
    let config = LotkaVolterraConfig::new(
        traditions,
        niche.overlap_matrix,
        0.01,
        SolverMethod::RK4,
    ).unwrap();

    let result = solve(&config, 5000).unwrap();
    let final_pops = result.trajectory.last().unwrap();
    let bio = biodiversity::compute(final_pops);

    assert!(bio.shannon_index >= 0.0);
    assert!(bio.simpson_index >= 0.0);
    assert!(bio.species_count >= 1);
    assert!(bio.evenness >= 0.0);
}

#[test]
fn single_tradition_logistic_growth() {
    let traditions = vec![make_tradition("Solo", 10.0, 1.0, 100.0, vec![0.0])];
    let config = LotkaVolterraConfig::new(
        traditions,
        vec![vec![1.0]],
        0.01,
        SolverMethod::RK4,
    ).unwrap();
    let result = solve(&config, 10000).unwrap();
    let final_pop = result.trajectory.last().unwrap()[0];
    assert!((final_pop - 100.0).abs() < 0.5, "solo tradition should reach K=100, got {final_pop}");
}

#[test]
fn three_way_coexistence() {
    let traditions = vec![
        make_tradition("Jazz", 30.0, 0.5, 100.0, vec![0.0, 0.0]),
        make_tradition("Blues", 30.0, 0.5, 100.0, vec![10.0, 0.0]),
        make_tradition("Rock", 30.0, 0.5, 100.0, vec![5.0, 8.66]),
    ];
    // Use low-competition matrix for coexistence
    let config = LotkaVolterraConfig::new(
        traditions,
        vec![vec![1.0, 0.1, 0.0], vec![0.1, 1.0, 0.1], vec![0.0, 0.1, 1.0]],
        0.01,
        SolverMethod::RK4,
    ).unwrap();
    let result = solve(&config, 50000).unwrap();
    let final_pops = result.trajectory.last().unwrap();
    for (i, &p) in final_pops.iter().enumerate() {
        assert!(p > 1.0, "tradition {i} went extinct: {p}");
    }
}
