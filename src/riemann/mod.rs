//! Riemann solver for the 1D Euler equations
//! following the method of Gottlieb and Groth.

pub mod waves;
use waves::{Contact, EulerSolution, Rarefaction, Shock};

use std::collections::VecDeque;

/// Euler equations state vector  
#[derive(Debug, Clone, Copy)]
pub struct EulerState {
    /// density in kg/m3
    pub density: f64,
    /// velocity in m/s
    pub velocity: f64,
    /// pressure in Pascals N/m2
    pub pressure: f64,
    /// Ratio of specific heats
    pub gamma: f64,
    /// Which side of the discontinuty the state is on
    pub side: StateSide,
}

/// The type of calculation
#[derive(Debug, Clone, Copy)]
pub enum StateSide {
    Left,
    Right,
}

/// The domain and interface between two initial states.
#[derive(Debug, Clone, Copy)]
pub struct DomainBounds {
    /// Initial interface point (x_0).
    pub interface: f64,
    /// Left-most point.
    pub left: f64,
    /// Right-most point.
    pub right: f64,
}

impl EulerState {
    /// Calculate the sound speed for a given state.
    pub fn sound_speed(&self) -> f64 {
        (self.gamma * self.pressure / self.density).sqrt()
    }

    /// Riemann invariant on the cell state.  
    pub fn big_gamma(&self) -> f64 {
        let sound_speed_term = self.sound_speed() * (2.0 / (self.gamma - 1.0));

        match self.side {
            StateSide::Left => self.velocity + sound_speed_term,
            StateSide::Right => self.velocity - sound_speed_term,
        }
    }
}

/// Exactly solve a Riemann problem between two Euler states iteratively.
pub fn solve_euler(left_ic: EulerState, right_ic: EulerState) -> EulerSolution {
    // check if there's a vacuum
    // first try was always true? probably written down wrong
    // if left_ic.big_gamma_left() - right_ic.big_gamma_right() > 0.0 {
    if left_ic.big_gamma() < 0.0 || right_ic.big_gamma() > 0.0 {
        return EulerSolution::RVR;
    }

    // initial guess
    let mut u_guess = velocity_guess(&left_ic, &right_ic);

    // stop criterion
    const EPSILON: f64 = 10e-6;

    // result struct
    let mut soln: EulerSolution;

    // pressure difference ratio that approaches 1.0
    let mut pressure_ratio: f64;

    // quick and dirty update polymorphism
    macro_rules! newton_update {
        ($left:ident, $right:ident, $guess:ident) => {
            $guess - ($left.pressure - $right.pressure) / ($left.dp_du - $right.dp_du)
        };
    }

    // carry out Newton's method
    loop {
        // check what the current solution is and update the guess accordingly
        let left_shock = u_guess < left_ic.velocity;
        let right_shock = u_guess > right_ic.velocity;

        soln = match (left_shock, right_shock) {
            // shock, shock
            (true, true) => {
                let left_state = Shock::create_shock(&left_ic, u_guess);
                let right_state = Shock::create_shock(&right_ic, u_guess);

                pressure_ratio = left_state.pressure / right_state.pressure;
                u_guess = newton_update!(left_state, right_state, u_guess);

                EulerSolution::SCS(left_state, Contact { velocity: u_guess }, right_state)
            }

            // rarefaction, shock
            (false, true) => {
                let left_state = Rarefaction::create_rarefaction(&left_ic, u_guess);
                let right_state = Shock::create_shock(&right_ic, u_guess);

                pressure_ratio = left_state.pressure / right_state.pressure;
                u_guess = newton_update!(left_state, right_state, u_guess);

                EulerSolution::RCS(left_state, Contact { velocity: u_guess }, right_state)
            }

            // shock, rarefaction
            (true, false) => {
                let left_state = Shock::create_shock(&left_ic, u_guess);
                let right_state = Rarefaction::create_rarefaction(&right_ic, u_guess);

                pressure_ratio = left_state.pressure / right_state.pressure;
                u_guess = newton_update!(left_state, right_state, u_guess);

                EulerSolution::SCR(left_state, Contact { velocity: u_guess }, right_state)
            }

            // rarefaction, rarefaction
            (false, false) => {
                let left_state = Rarefaction::create_rarefaction(&left_ic, u_guess);
                let right_state = Rarefaction::create_rarefaction(&right_ic, u_guess);

                pressure_ratio = left_state.pressure / right_state.pressure;
                u_guess = newton_update!(left_state, right_state, u_guess);

                EulerSolution::RCR(left_state, Contact { velocity: u_guess }, right_state)
            }
        };

        // check if we can stop
        if (1.0 - pressure_ratio).abs() < EPSILON {
            break soln; // return whichever solution it was
        }
    }
}

/// Give an initial guess for the velocity following Groth and Gottlieb.
fn velocity_guess(left: &EulerState, right: &EulerState) -> f64 {
    (left.big_gamma() * z_term(left, right) + right.big_gamma()) / (1.0 + z_term(left, right))
}

/// The z-term depending on the left and right states.
fn z_term(left: &EulerState, right: &EulerState) -> f64 {
    // find sigma term
    let sigma = if left.pressure <= right.pressure {
        left.gamma
    } else {
        right.gamma
    };

    // calculate z
    let exponent = (sigma - 1.0) / (2.0 * sigma);

    ((left.gamma - 1.0) / (right.gamma - 1.0)
        * (right.sound_speed() / left.sound_speed())
        * (left.pressure / right.pressure))
        .powf(exponent)
}

#[cfg(test)]
mod tests {

    use super::*;

    extern crate approx;
    use approx::assert_relative_eq;

    #[test]
    fn test_sound_speed() {
        let left_state = EulerState {
            density: 2.281,
            velocity: 164.83,
            pressure: 201.17e3,
            gamma: 1.4,
        };

        assert_relative_eq!(left_state.sound_speed(), 351.3848, epsilon = 0.001);
    }
}
