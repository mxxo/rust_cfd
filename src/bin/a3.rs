//!  Finite-Volume methods assignment 3
//!  Max Orok, November-December 2019
//! -----------------------------------------------------------------------------
//!   Approximate flux function comparison for the 1D Euler equations
//! -----------------------------------------------------------------------------

use rust_cfd::euler::{EulerSolution1d, PrimitiveResult, PrimitiveState};

// assignment 3
use rust_cfd::fluxes;
use rust_cfd::fluxes::FluxFunction;

use rust_cfd::riemann::DomainBounds;

// given as the same for both gases
const GAMMA: f64 = 1.4;

struct SolnSpec {
    pub t_final: f64,
    pub bounds: DomainBounds,
    pub num_cells: usize,
    pub cfl: f64,
}

fn main() {

    let flux_fn = fluxes::Exact {};
//     let flux_fn = fluxes::Roe {};
//     let flux_fn = fluxes::RoeEntropyFix {};
//     let flux_fn = fluxes::Hlle {};

    let num_cells = 1000;
    let cfl = 0.99;

    println!("\ncase 1\n");
    case_1(num_cells, cfl, flux_fn);
    println!("\ncase 2\n");
    case_2(num_cells, cfl, flux_fn);
    println!("\ncase 3\n");
    case_3(num_cells, cfl, flux_fn);
    println!("\ncase 4\n");
    case_4(num_cells, cfl, flux_fn);
}

fn case_1(num_cells: usize, cfl: f64, flux_fn: impl FluxFunction) {
    let left_state = PrimitiveState {
        density: 2.281,
        velocity: 164.83,
        pressure: 201.17e3,
        gamma: GAMMA,
    };

    let right_state = PrimitiveState {
        density: 1.408,
        velocity: 0.0,
        pressure: 101.1e3,
        gamma: GAMMA,
    };

    let t_final = 12e-3;

    let bounds = DomainBounds {
        left: 0.0,
        interface: 2.0,
        right: 10.0,
    };

    let soln_spec = SolnSpec {
        t_final,
        bounds,
        num_cells,
        cfl,
    };

    solve_euler_eqns(left_state, right_state, flux_fn, soln_spec);
}

fn case_2(num_cells: usize, cfl: f64, flux_fn: impl FluxFunction) {
    let left_state = PrimitiveState {
        density: 1.045,
        velocity: 200.0,
        pressure: 300e3,
        gamma: GAMMA,
    };

    let right_state = PrimitiveState {
        density: 3.483,
        velocity: 200.0,
        pressure: 300e3,
        gamma: GAMMA,
    };

    let t_final = 25e-3;

    let bounds = DomainBounds {
        left: 0.0,
        interface: 2.0,
        right: 10.0,
    };

    let soln_spec = SolnSpec {
        t_final,
        bounds,
        num_cells,
        cfl,
    };

    solve_euler_eqns(left_state, right_state, flux_fn, soln_spec);
}

fn case_3(num_cells: usize, cfl: f64, flux_fn: impl FluxFunction) {
    let left_state = PrimitiveState {
        density: 1.598,
        velocity: -383.64,
        pressure: 91.88e3,
        gamma: GAMMA,
    };

    let right_state = PrimitiveState {
        density: 2.787,
        velocity: -216.97,
        pressure: 200e3,
        gamma: GAMMA,
    };

    let t_final = 35e-3; // s

    let bounds = DomainBounds {
        left: 0.0,
        interface: 5.0,
        right: 10.0,
    };

    let soln_spec = SolnSpec {
        t_final,
        bounds,
        num_cells,
        cfl,
    };

    solve_euler_eqns(left_state, right_state, flux_fn, soln_spec);
}

fn case_4(num_cells: usize, cfl: f64, flux_fn: impl FluxFunction) {
    let left_state = PrimitiveState {
        density: 4.696,
        velocity: 0.0,
        pressure: 404.4e3,
        gamma: GAMMA,
    };

    let right_state = PrimitiveState {
        density: 1.408,
        velocity: 0.0,
        pressure: 101.1e3,
        gamma: GAMMA,
    };

    let t_final = 7e-3; // s
    let bounds = DomainBounds {
        left: 0.0,
        interface: 5.0,
        right: 10.0,
    };

    let soln_spec = SolnSpec {
        t_final,
        bounds,
        num_cells,
        cfl,
    };

    solve_euler_eqns(left_state, right_state, flux_fn, soln_spec);
}

fn solve_euler_eqns(
    left_state: PrimitiveState,
    right_state: PrimitiveState,
    flux_fn: impl FluxFunction,
    soln_spec: SolnSpec,
) {
    let soln = EulerSolution1d::init(
        left_state,
        right_state,
        soln_spec.bounds,
        soln_spec.num_cells,
        soln_spec.cfl,
    );

    // let b_conds = FixedBoundary {
    //     left_state,
    //     right_state,
    // };

    let res = soln.first_order_time_march(flux_fn, soln_spec.t_final);

    let res: Vec<_> = res
        .into_iter()
        .map(|cell| cell.to_primitive_result())
        .collect();

    println!("coord\tdensity\tvelocity\tpressure");

    for PrimitiveResult { state, coord, .. } in res {
        println!(
            "{:.4}\t{:.4}\t{:.4}\t{:.4}",
            coord, state.density, state.velocity, state.pressure
        );
    }
}
