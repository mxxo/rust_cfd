//!  Finite-Volume methods assignment 2
//!  Max Orok, October-November 2019
//! -----------------------------------------------------------------------------
//!   Iterative Riemann solver for the 1D Euler equations
//! -----------------------------------------------------------------------------

use rust_cfd::riemann::{solve_euler, EulerState, StateSide};

fn main() {
    // given as the same for both sides
    const GAMMA: f64 = 1.4;

    // case 1
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
    soln.reconstruct(t_final);

    // case 2
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
    soln.reconstruct(t_final);

    // case 3
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
    soln.reconstruct(t_final);

    // case 4
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
    soln.reconstruct(t_final);
}
