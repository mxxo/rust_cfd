//!  Finite-Volume methods assignment 3
//!  Max Orok, November-December 2019
//! -----------------------------------------------------------------------------
//!   Approximate flux function comparison for the 1D Euler equations
//! -----------------------------------------------------------------------------

use rust_cfd::riemann::waves::DataPoint;
use rust_cfd::riemann::{solve_euler, DomainBounds, EulerState, StateSide};

use rust_cfd::euler::{FluxType, EulerSolution1d};
//
// given as the same for both gases
const GAMMA: f64 = 1.4;

fn main() {

    let flux_type = FluxType::Exact;
    let num_cells = 100;

    case_1(flux_type, num_cells);
    case_2(flux_type, num_cells);
    case_3(flux_type, num_cells);
    case_4(flux_type, num_cells);

}

fn case_1(flux: FluxType, num_cells: u32) {

    let left_state = EulerState {
        density: 2.281,
        velocity: 164.83,
        pressure: 201.17e3,
        gamma: GAMMA,
        side: StateSide::Left,
    };

    let right_state = EulerState {
        density: 1.408,
        velocity: 0.0,
        pressure: 101.1e3,
        gamma: GAMMA,
        side: StateSide::Right,
    };

    let soln = solve_euler(left_state, right_state);
    dbg!(&soln);

    let t_final = 12e-3; // s
    let vec_out = soln.reconstruct(
        left_state,
        right_state,
        DomainBounds {
            left: 0.0,
            interface: 2.0,
            right: 10.0,
        },
        t_final,
    );

    println!("<< case 1 >>\n");
    println!("coord\tdensity\tvelocity\tpressure");

    for DataPoint {
        coord,
        density,
        velocity,
        pressure,
    } in vec_out
    {
        println!(
            "{:.4}\t{:.4}\t{:.4}\t{:.4}",
            coord, density, velocity, pressure
        );
    }
}

fn case_2(flux: FluxType, num_cells: u32) {
    let left_state = EulerState {
        density: 1.045,
        velocity: 200.0,
        pressure: 300e3,
        gamma: GAMMA,
        side: StateSide::Left,
    };
    let right_state = EulerState {
        density: 3.483,
        velocity: 200.0,
        pressure: 300e3,
        gamma: GAMMA,
        side: StateSide::Right,
    };

    let soln = solve_euler(left_state, right_state);
    dbg!(&soln);

    let t_final = 25e-3; // s
    let vec_out = soln.reconstruct(
        left_state,
        right_state,
        DomainBounds {
            left: 0.0,
            interface: 2.0,
            right: 10.0,
        },
        t_final,
    );

    println!("<< case 2 >>\n");
    println!("coord\tdensity\tvelocity\tpressure");

    for DataPoint {
        coord,
        density,
        velocity,
        pressure,
    } in vec_out
    {
        println!(
            "{:.4}\t{:.4}\t{:.4}\t{:.4}",
            coord, density, velocity, pressure
        );
    }
}

fn case_3(flux: FluxType, num_cells: u32) {
    let left_state = EulerState {
        density: 1.598,
        velocity: -383.64,
        pressure: 91.88e3,
        gamma: GAMMA,
        side: StateSide::Left,
    };
    let right_state = EulerState {
        density: 2.787,
        velocity: -216.97,
        pressure: 200e3,
        gamma: GAMMA,
        side: StateSide::Right,
    };

    let soln = solve_euler(left_state, right_state);
    dbg!(&soln);

    let t_final = 35e-3; // s
    let vec_out = soln.reconstruct(
        left_state,
        right_state,
        DomainBounds {
            left: 0.0,
            interface: 5.0,
            right: 10.0,
        },
        t_final,
    );
    println!("<< case 3 >>\n");
    println!("coord\tdensity\tvelocity\tpressure");

    for DataPoint {
        coord,
        density,
        velocity,
        pressure,
    } in vec_out
    {
        println!(
            "{:.4}\t{:.4}\t{:.4}\t{:.4}",
            coord, density, velocity, pressure
        );
    }
}

fn case_4(flux: FluxType, num_cells: u32) {
    let left_state = EulerState {
        density: 4.696,
        velocity: 0.0,
        pressure: 404.4e3,
        gamma: GAMMA,
        side: StateSide::Left,
    };
    let right_state = EulerState {
        density: 1.408,
        velocity: 0.0,
        pressure: 101.1e3,
        gamma: GAMMA,
        side: StateSide::Right,
    };

    let soln = solve_euler(left_state, right_state);
    dbg!(&soln);

    let t_final = 7e-3; // s
    let vec_out = soln.reconstruct(
        left_state,
        right_state,
        DomainBounds {
            left: 0.0,
            interface: 5.0,
            right: 10.0,
        },
        t_final,
    );
    println!("<< case 4 >>\n");
    println!("coord\tdensity\tvelocity\tpressure");

    for DataPoint {
        coord,
        density,
        velocity,
        pressure,
    } in vec_out
    {
        println!(
            "{:.4}\t{:.4}\t{:.4}\t{:.4}",
            coord, density, velocity, pressure
        );
    }
}
