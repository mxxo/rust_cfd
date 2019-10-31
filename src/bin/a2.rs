//!  Finite-Volume methods assignment 2
//!  Max Orok, October 2019
//! -----------------------------------------------------------------------------
//! Iterative Riemann solver
//! -----------------------------------------------------------------------------
//!  There are five possible wave patterns for 1D Euler
//!
//!  S - Shock
//!  C - Contact surface
//!  R - Rarefaction/Expansion wave
//!  V - Vacuum
//!
//!  1. SCS
//!  2. RCS
//!  3. SCR
//!  4. RCR
//!  5. RVR (very rare)
//!
use rust_cfd::riemann::{solve_euler, EulerState};

fn main() {
    
    const GAMMA: f64 = 1.4; 

    // case 1
    let left_state = EulerState {
        density: 2.281,
        velocity: 164.83,
        pressure: 201.17e3,
        gamma: GAMMA
    };

    let right_state = EulerState {
        density: 1.408,
        velocity: 0.0,
        pressure: 101.1e3,
        gamma: GAMMA
    };
    
    let t_final = 12e-3; // s 

    dbg!(&left_state);
    dbg!(&right_state);
    
    solve_euler(left_state, right_state, t_final); 

    // case 2
    let left_state = EulerState {
        density: 1.045,
        velocity: 200.0,
        pressure: 300e3,
        gamma: GAMMA, 
    };
    let right_state = EulerState {
        density: 3.483,
        velocity: 200.0,
        pressure: 300e3,
        gamma: GAMMA, 
    };

    let t_final = 25e-3; // s

    dbg!(&left_state);
    dbg!(&right_state);

    solve_euler(left_state, right_state, t_final); 
    
    // case 3
    let left_state = EulerState {
        density: 1.598,
        velocity: -383.64,
        pressure: 91.88e3,
        gamma: GAMMA, 
    };
    let right_state = EulerState {
        density: 2.787,
        velocity: -216.97,
        pressure: 200e3,
        gamma: GAMMA, 
    };

    let t_final = 35e-3; // s

    dbg!(&left_state);
    dbg!(&right_state);

    solve_euler(left_state, right_state, t_final); 

    // case 4
    let left_state = EulerState {
        density: 4.696,
        velocity: 0.0,
        pressure: 404.4e3,
        gamma: GAMMA, 
    };
    let right_state = EulerState {
        density: 1.408,
        velocity: 0.0,
        pressure: 101.1e3,
        gamma: GAMMA, 
    };

    let t_final = 7e-3; // s

    dbg!(&left_state);
    dbg!(&right_state);

    solve_euler(left_state, right_state, t_final); 
}
