//! Riemann solver for the 1D Euler equations
//! following the method of Gottlieb and Groth.

extern crate approx;
use approx::assert_relative_eq;

const EPSILON: f64 = 10e-6;

/// Euler equations state vector  
#[derive(Debug)]
pub struct EulerState {
    /// density in kg/m3
    pub density: f64,
    /// velocity in m/s
    pub velocity: f64,
    /// pressure in Pascals N/m2
    pub pressure: f64,
    /// Ratio of specific heats
    pub gamma: f64,
}

impl EulerState {
    /// Calculate the sound speed for a given state.
    pub fn sound_speed(&self) -> f64 {
        (self.gamma * self.pressure / self.density).sqrt()
    }

    /// Riemann invariant on the left state.  
    pub fn big_gamma_left(&self) -> f64 {
        self.velocity + self.sound_speed() * (2.0 / (self.gamma - 1.0))
    }

    /// Riemann invariant on the right state.
    pub fn big_gamma_right(&self) -> f64 {
        self.velocity - self.sound_speed() * (2.0 / (self.gamma - 1.0))
    }
}

/// Solve the Euler equations iteratively between the left and right states.
pub fn solve_euler(left_ic: EulerState, right_ic: EulerState, t_final: f64) {
    // check if there's a vacuum
    if left_ic.big_gamma_left() - right_ic.big_gamma_right() > 0.0 {
        panic!("Vacuum detected, aborting");
    }

    // initial guess
    let u_0 = velocity_guess(&left_ic, &right_ic);

    // begin Newton's method
}

/// Async Newton's method for 1D Euler equations (some sort of async Stream?)
// select the first successful future?
// https://rust-lang-nursery.github.io/futures-api-docs/0.3.0-alpha.19/futures/future/fn.select_ok.html

/// Give an initial guess for the velocity following Groth and Gottlieb.
fn velocity_guess(left: &EulerState, right: &EulerState) -> f64 {
    (left.big_gamma_left() * z_term(left, right) + right.big_gamma_right())
        / (1.0 + z_term(left, right))
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
