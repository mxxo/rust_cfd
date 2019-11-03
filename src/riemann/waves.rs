//! Waves and wave patterns for the Euler equations in 1D.

use super::*;

/// Euler solution wave patterns in 1D.  
///
/// There are five possible outcomes.
/// In practice, vacuum is extremely rare, and usually neglected.
///
/// The components are abbreviated as:
/// * S - Shock
/// * C - Contact surface
/// * R - Rarefaction/Expansion wave
/// * V - Vacuum
#[derive(Debug)]
pub enum EulerSolution {
    /// Shock, contact surface, shock
    SCS(Shock, Contact, Shock),
    /// Rarefaction, contact surface, shock
    RCS(Rarefaction, Contact, Shock),
    /// Shock, contact surface,
    SCR(Shock, Contact, Rarefaction),
    /// Rarefaction, contact surface, rarefaction
    RCR(Rarefaction, Contact, Rarefaction),
    /// Rarefaction, vaccum, rarefaction
    RVR,
}

/// A contact surface between two gases.
#[derive(Debug)]
pub struct Contact {
    /// Contact surface velocity (m/s).
    pub velocity: f64,
}

/// A shock wave in a gas.  
#[derive(Debug)]
pub struct Shock {
    /// The shock mach number (dimensionless).
    pub mach_number: f64,
    /// Pressure in Pascals (Pa).
    pub pressure: f64,
    /// Pressure derivative relative to the contact surface velocity.
    pub dp_du: f64,
}

impl Shock {
    /// Create a shock for a given state.
    pub fn create_shock(state: &EulerState, velocity_guess: f64) -> Self {
        let mach_number = Self::mach(state, velocity_guess);
        let pressure = Self::pressure(state, velocity_guess, mach_number);
        let dp_du = Self::pressure_derivative(state, mach_number);

        Self {
            mach_number,
            pressure,
            dp_du,
        }
    }

    /// Shock Mach number for a given velocity.
    fn mach(state: &EulerState, velocity_guess: f64) -> f64 {
        let gamma_term =
            ((state.gamma + 1.0) / 4.0) * ((velocity_guess - state.velocity) / state.sound_speed());

        let sqrt_term = (1.0 + gamma_term * gamma_term).sqrt();

        match state.side {
            StateSide::Left => gamma_term - sqrt_term,
            StateSide::Right => gamma_term + sqrt_term,
        }
    }

    /// Shock pressure for a given velocity.
    fn pressure(state: &EulerState, velocity_guess: f64, mach_number: f64) -> f64 {
        let mach_term =
            state.gamma / state.sound_speed() * (velocity_guess - state.velocity) * mach_number;
        state.pressure * (1.0 + mach_term)
    }

    /// Shock pressure derivative relative to velocity for a given velocity.
    fn pressure_derivative(state: &EulerState, mach_number: f64) -> f64 {
        let mach_cubed = mach_number * mach_number * mach_number;
        let mach_denom = 1.0 + mach_number * mach_number;

        2.0 * state.gamma * state.pressure / state.sound_speed() * mach_cubed / mach_denom
    }
}

/// A rarefaction wave in a gas.  
#[derive(Debug)]
pub struct Rarefaction {
    /// Sound speed for the rarefaction wave (m/s).
    pub sound_speed: f64,
    /// Pressure in Pascals (Pa).
    pub pressure: f64,
    /// Pressure derivative relative to the contact surface velocity.
    pub dp_du: f64,
}

impl Rarefaction {
    /// Create a rarefaction wave for a given state.
    pub fn create_rarefaction(state: &EulerState, velocity_guess: f64) -> Self {
        let sound_speed = Self::sound_speed(state, velocity_guess);
        let pressure = Self::pressure(state, sound_speed);
        let dp_du = Self::pressure_derivative(state, sound_speed, pressure);

        Self {
            sound_speed,
            pressure,
            dp_du,
        }
    }

    /// Sound speed in the rarefaction wave.
    fn sound_speed(state: &EulerState, velocity_guess: f64) -> f64 {
        let gamma_term =
            (state.gamma - 1.0) / 2.0 * (velocity_guess - state.velocity) / state.sound_speed();

        match state.side {
            StateSide::Left => state.sound_speed() * (1.0 - gamma_term),
            StateSide::Right => state.sound_speed() * (1.0 + gamma_term),
        }
    }

    /// Rarefaction pressure for a given velocity.
    fn pressure(state: &EulerState, sound_speed: f64) -> f64 {
        let gamma_expo = 2.0 * state.gamma / (state.gamma - 1.0);
        state.pressure * (sound_speed / state.sound_speed()).powf(gamma_expo)
    }

    /// Rarefaction pressure derivative relative to velocity for a given velocity.
    fn pressure_derivative(state: &EulerState, sound_speed: f64, pressure: f64) -> f64 {
        match state.side {
            StateSide::Left => -state.gamma * pressure / sound_speed,
            StateSide::Right => state.gamma * pressure / sound_speed,
        }
    }
}
