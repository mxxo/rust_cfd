//! Riemann solver for the 1D Euler equations
//! following the method of Gottlieb and Groth.

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
    
    /// Calculate the sound speed for a given state 
    pub fn sound_speed(&self) -> f64 {
        (self.gamma * self.pressure / self.density).sqrt()
    }

}

/// Solve the Euler equations iteratively between the left and right states
pub fn solve_euler(left: EulerState, right: EulerState, t_final: f64) { 

}

// /// Sound speed equation 
// pub fn sound_speed(pressure: f64, density: f64, gamma: f64) -> f64 { 
//     (gamma * pressure / density).sqrt()
// }

/// Riemann invariant on the left state  
fn big_gamma_left(velocity: f64, gamma: f64, sound_speed: f64) -> f64 {
    velocity + sound_speed * (2.0 / (gamma - 1.0))
}

/// Riemann invariant on the right state  
fn big_gamma_right(velocity: f64, gamma: f64, sound_speed: f64) -> f64 {
    velocity - sound_speed * (2.0 / (gamma - 1.0))
}
