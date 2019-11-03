//! Riemann solver for the 1D Euler equations
//! following the method of Gottlieb and Groth.

extern crate approx;
use approx::assert_relative_eq;

pub mod waves;
use waves::{Contact, EulerSolution, Rarefaction, Shock};

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
    /// Which side of the discontinuty the state is on
    pub side: StateSide,
}

/// The type of calculation
#[derive(Debug)]
pub enum StateSide {
    Left,
    Right,
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

/// Solve a Riemann problem between two Euler states iteratively.
pub fn solve_euler(left_ic: EulerState, right_ic: EulerState, t_final: f64) -> EulerSolution {
    dbg!(left_ic.big_gamma());
    dbg!(right_ic.big_gamma());

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

    let mut soln: EulerSolution;

    // pressure difference ratio that approaches 1.0
    let mut pressure_ratio: f64;

    // carry out Newton's method
    loop {
        // check what the current solution is and update the guess accordingly
        soln = match (u_guess < left_ic.velocity, u_guess > right_ic.velocity) {
            
            // shock, shock
            (true, true) => {
                let left_state = Shock::create_shock(&left_ic, u_guess);
                let right_state = Shock::create_shock(&right_ic, u_guess);
                pressure_ratio = left_state.pressure / right_state.pressure;

                u_guess = newton_update(
                    u_guess,
                    left_state.pressure,
                    right_state.pressure,
                    left_state.dp_du,
                    right_state.dp_du,
                );

                EulerSolution::SCS(left_state, Contact { velocity: u_guess }, right_state)
            }
            
            // rarefaction, shock
            (false, true) => { 
                let left_state = Rarefaction::create_rarefaction(&left_ic, u_guess);
                let right_state = Shock::create_shock(&right_ic, u_guess);
                pressure_ratio = left_state.pressure / right_state.pressure;

                u_guess = newton_update(
                    u_guess,
                    left_state.pressure,
                    right_state.pressure,
                    left_state.dp_du,
                    right_state.dp_du,
                );

                EulerSolution::RCS(left_state, Contact { velocity: u_guess }, right_state)

            },
            
            // shock, rarefaction
            (true, false) => { 
                let left_state = Shock::create_shock(&left_ic, u_guess);
                let right_state = Rarefaction::create_rarefaction(&right_ic, u_guess);
                pressure_ratio = left_state.pressure / right_state.pressure;

                u_guess = newton_update(
                    u_guess,
                    left_state.pressure,
                    right_state.pressure,
                    left_state.dp_du,
                    right_state.dp_du,
                );

                EulerSolution::SCR(left_state, Contact { velocity: u_guess }, right_state)
            },
            
            // rarefaction, rarefaction
            (false, false) => { 
                let left_state = Rarefaction::create_rarefaction(&left_ic, u_guess);
                let right_state = Rarefaction::create_rarefaction(&right_ic, u_guess);
                pressure_ratio = left_state.pressure / right_state.pressure;

                u_guess = newton_update(
                    u_guess,
                    left_state.pressure,
                    right_state.pressure,
                    left_state.dp_du,
                    right_state.dp_du,
                );

                EulerSolution::RCR(left_state, Contact { velocity: u_guess }, right_state)
            }, 
        };

        // check if we can stop
        if (1.0 - pressure_ratio).abs() < EPSILON {
            break soln; // return whichever solution it was 
        }
    }
}

/// Newton's method to solve Euler states
fn newton_update(
    current_guess: f64,
    left_pressure: f64,
    right_pressure: f64,
    left_dp_du: f64,
    right_dp_du: f64,
) -> f64 {
    current_guess - (left_pressure - right_pressure) / (left_dp_du - right_dp_du)
}

/// Async Newton's method for 1D Euler equations (some sort of async Stream?)
// select the first successful future?
// https://rust-lang-nursery.github.io/futures-api-docs/0.3.0-alpha.19/futures/future/fn.select_ok.html

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
