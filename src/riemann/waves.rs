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
    ) -> VecDeque<DataPoint> {
        // find the wave pattern solution around an interface at x = 0
        // comes back sorted
        let mut soln = self.extrapolate(left, right, time);

        // shift everyone over by the correct amount
        for data_point in &mut soln {
            data_point.coord += bounds.interface;
        }

        // switch to deque for easy inserts at the front
        let mut soln = VecDeque::from(soln);

        // check whether left state is still on the domain and
        // and if so add points at the start using the left state
        if bounds.left < soln.front().unwrap().coord {
            soln.push_front(DataPoint {
                coord: bounds.left,
                value: left.density,
            });
        }

        // check whether right state is still on the domain and
        // and if so add points at the end using the right state
        if bounds.right > soln.back().unwrap().coord {
            soln.push_back(DataPoint {
                coord: bounds.right,
                value: right.density,
            });
        }

        soln
    }

    /// Extrapolate the solution and return a list of data points sorted by coordinate.
    fn extrapolate(&self, left: EulerState, right: EulerState, time: f64) -> Vec<DataPoint> {
        // quick and dirty plot polymorphism
        macro_rules! plot_waves {
            ($left_wave:ident, $contact:ident, $right_wave:ident) => {{
                let mut res = Vec::new();

                // left wave data
                let left_data = $left_wave.plot(&left, $contact, time);
                // get state at rightmost point
                let left_state_star = *left_data.last().unwrap();
                res.extend(left_data);

                // right wave data
                let right_data = $right_wave.plot(&right, $contact, time);
                // get state at rightmost point
                let right_state_star = *right_data.last().unwrap();
                res.extend(right_data);

                // plot contact surface using left and right wave information
                res.extend($contact.plot(time, left_state_star, right_state_star));

                res.sort_by(|a, b| a.coord.partial_cmp(&b.coord).unwrap());
                res
            }};
        }

        match self {
            EulerSolution::RVR => {
                panic!("Unable to plot Rarefaction, Vacuum, Rarefaction solution")
            }
            EulerSolution::RCR(left_raref, contact, right_raref) => {
                plot_waves!(left_raref, contact, right_raref)
            }
            EulerSolution::RCS(left_raref, contact, right_shock) => {
                plot_waves!(left_raref, contact, right_shock)
            }
            EulerSolution::SCR(left_shock, contact, right_raref) => {
                plot_waves!(left_shock, contact, right_raref)
            }
            EulerSolution::SCS(left_shock, contact, right_shock) => {
                plot_waves!(left_shock, contact, right_shock)
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
/// A 1D data point.
pub struct DataPoint {
    /// Coordinate
    coord: f64,
    /// Value
    value: f64,
}

/// A contact surface between two gases.
#[derive(Debug, Clone, Copy)]
pub struct Contact {
    /// Contact surface velocity (m/s).
    pub velocity: f64,
}

impl Contact {
    pub fn plot(&self, time: f64, left_state: DataPoint, right_state: DataPoint) -> Vec<DataPoint> {
        // rough plot for now
        // put our known points in the result vector
        vec![
            DataPoint {
                coord: left_state.coord,
                value: left_state.value,
            },
            DataPoint {
                coord: self.extrapolate(time),
                value: left_state.value,
            },
            DataPoint {
                coord: right_state.coord,
                value: right_state.value,
            },
        ]
    }

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

    /// Plot shock data.
    pub fn plot(&self, state: &EulerState, contact: &Contact, time: f64) -> Vec<DataPoint> {
        let shock_coord = self.extrapolate(state, time);

        let shift = 0.0001;

        match state.side {
            StateSide::Left => {
                vec![
                    // left state data first
                    DataPoint {
                        coord: shock_coord,
                        value: state.density,
                    },
                    // contact surface data
                    DataPoint {
                        coord: shock_coord + shift,
                        value: Rarefaction::density(state, contact.velocity),
                    },
                ]
            }
            StateSide::Right => {
                vec![
                    // contact surface data
                    DataPoint {
                        coord: shock_coord,
                        value: Rarefaction::density(state, contact.velocity),
                    },
                    // right state data
                    DataPoint {
                        coord: shock_coord + shift,
                        value: Rarefaction::density(state, state.velocity),
                    },
                ]
            }
        }
    }

    /// Extrapolate a shock given an initial state.
    pub fn extrapolate(&self, state: &EulerState, time: f64) -> f64 {
        match state.side {
            StateSide::Left => time * (state.velocity - self.mach_number * state.sound_speed()),
            StateSide::Right => time * (state.velocity + self.mach_number * state.sound_speed()),
        }
    }

    // --
    // Shock relations
    // --

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

    /// Plot the rarefaction wave
    pub fn plot(&self, state: &EulerState, contact: &Contact, time: f64) -> Vec<DataPoint> {
        let (head_pos, tail_pos) = self.extrapolate(state, contact, time);

        // find head rarefaction wave's density given state's velocity
        let head_data = DataPoint {
            coord: head_pos,
            value: Self::density(state, state.velocity),
        };

        // find rarefaction wave tail's density given the contact surface
        let tail_data = DataPoint {
            coord: tail_pos,
            value: Self::density(state, contact.velocity),
        };

        if head_data.coord < tail_data.coord {
            vec![head_data, tail_data]
        } else {
            vec![tail_data, head_data]
        }
    }

    /// Extrapolate a rarefaction wave given the initial state and a contact surface.
    /// Returns the positions of the rarefaction head and tail as `(x_head, x_tail)`
    pub fn extrapolate(&self, state: &EulerState, contact: &Contact, time: f64) -> (f64, f64) {
        match state.side {
            StateSide::Left => {
                let head_pos = (state.velocity - state.sound_speed()) * time;
                let tail_pos = (contact.velocity - self.sound_speed) * time;
                (head_pos, tail_pos)
            }
            StateSide::Right => {
                let head_pos = state.velocity + state.sound_speed() * time;
                let tail_pos = contact.velocity + self.sound_speed * time;
                (head_pos, tail_pos)
            }
        }
    }

    // --
    // Rarefaction relations
    // --

    /// Rarefaction density for a given velocity.
    pub fn density(state: &EulerState, velocity: f64) -> f64 {
        let local_sound_speed = Self::sound_speed(state, velocity);
        let local_pressure = Self::pressure(state, local_sound_speed);

        state.gamma * local_pressure / (local_sound_speed * local_sound_speed)
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
