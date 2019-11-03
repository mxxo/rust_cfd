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
#[derive(Debug, Clone, Copy)]
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

impl EulerSolution {
    /// Plot the exact solution at a given time.  
    pub fn reconstruct(
        &self,
        left: EulerState,
        right: EulerState,
        bounds: DomainBounds,
        time: f64,
    ) -> Vec<(f64, f64)> {
        match self {
            EulerSolution::RVR => {
                panic!("Unable to plot Rarefaction, Vacuum, Rarefaction solution")
            }
            EulerSolution::RCR(left_raref, contact, right_raref)  => { 
                // left rarefaction
                let (l_head, l_tail) = left_raref.extrapolate(&left, contact, time); 
                // contact 
                let contact_pos = contact.extrapolate(time); 
                // right rarefaction
                let (r_head, r_tail) = right_raref.extrapolate(&right, contact, time); 

                vec![(l_head, 0.0), (l_tail, 0.0), (contact_pos, 0.0), (r_tail, 0.0), (r_head, 0.0)]
            },
            EulerSolution::RCS(left_raref, contact, right_shock)  => { 
                // left rarefaction
                let (l_head, l_tail) = left_raref.extrapolate(&left, contact, time); 
                // contact 
                let contact_pos = contact.extrapolate(time); 
                // right shock
                let shock_pos = right_shock.extrapolate(&right, time); 

                vec![(l_head, 0.0), (l_tail, 0.0), (contact_pos, 0.0), (shock_pos, 0.0)]
            },
            
            // some sort of tagging system for what type of data this is? 

            _ => unimplemented!(), 
        }
    }
}

/// A contact surface between two gases.
#[derive(Debug, Clone, Copy)]
pub struct Contact {
    /// Contact surface velocity (m/s).
    pub velocity: f64,
}

impl Contact {
    pub fn extrapolate(&self, time: f64) -> f64 {
        self.velocity * time
    }
}

/// A shock wave in a gas.  
#[derive(Debug, Clone, Copy)]
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

    /// Extrapolate a shock given an initial state.
    pub fn extrapolate(&self, state: &EulerState, time: f64) -> f64 {
        match state.side {
            StateSide::Left => time * (state.velocity - self.mach_number * state.sound_speed()),
            StateSide::Right => time * (state.velocity + self.mach_number * state.sound_speed()),
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
#[derive(Debug, Clone, Copy)]
pub struct Rarefaction {
    /// Sound speed for the rarefaction wave (m/s).
    pub sound_speed: f64,
    /// Pressure in Pascals (Pa).
    pub pressure: f64,
    /// Pressure derivative relative to the contact surface velocity.
    pub dp_du: f64,
}

impl Rarefaction {

    /// Rarefaction density for a given velocity. 
    pub fn density(&self, state: &EulerState, velocity: f64) -> f64 {
        let local_sound_speed = Self::sound_speed(state, velocity); 
        let local_pressure = Self::pressure(state, local_sound_speed); 
        
        state.gamma * local_pressure / (local_sound_speed * local_sound_speed)
    }

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
    
    /// Extrapolate a rarefaction wave given the initial state and a contact surface.
    /// Returns the positions of the rarefaction head and tail as `(x_head, x_tail)`
    pub fn extrapolate(&self, state: &EulerState, contact: &Contact, time: f64) -> (f64, f64) { 
        match state.side {
            StateSide::Left => { 
                let head = state.velocity - state.sound_speed(); 
                let tail = contact.velocity - self.sound_speed; 
                (head * time, tail * time)
            },
            StateSide::Right => { 
                let head = state.velocity + state.sound_speed(); 
                let tail = contact.velocity + self.sound_speed; 
                (head * time, tail * time)
            }, 
        }
    }
    
    // -- 
    // Rarefaction relations 
    // -- 

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
