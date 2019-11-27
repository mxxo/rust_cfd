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

        // only keep those points inside the bounds
        soln.retain(|&data_point| {
            data_point.coord > bounds.left && data_point.coord < bounds.right
        });

        // switch to deque for easy inserts at the front
        let mut soln = VecDeque::from(soln);

        // -- tag on start and end points

        // check whether left state is still on the domain and
        // and if so add points at the start using the left state
        if bounds.left < soln.front().unwrap().coord {
            soln.push_front(DataPoint {
                coord: bounds.left,
                density: left.density,
                velocity: left.velocity,
                pressure: left.pressure,
            });
        }

        // check whether right state is still on the domain and
        // and if so add points at the end using the right state
        if bounds.right > soln.back().unwrap().coord {
            soln.push_back(DataPoint {
                coord: bounds.right,
                density: right.density,
                velocity: right.velocity,
                pressure: right.pressure,
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
                // get state at rightmost point -- will be the left_star's data
                let left_state_star = *left_data.last().unwrap();
                res.extend(left_data);


                // right wave data
                let right_data = $right_wave.plot(&right, $contact, time);
                // get state at leftmost point -- the right_star's data
                let right_state_star = right_data[0];

                res.extend(right_data);

//                dbg!(&res);


                // plot contact surface using left and right wave information
                res.extend($contact.plot(time, left_state_star, right_state_star));

 //               dbg!(&res);

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
/// A 1D primitive state variables data point.
pub struct DataPoint {
    /// Coordinate
    pub coord: f64,
    /// Density
    pub density: f64,
    /// Velocity
    pub velocity: f64,
    /// Pressure
    pub pressure: f64,
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

        let shift = 0.0001;

        // put our known points in the result vector
        vec![
            DataPoint {
                coord: left_state.coord,
                density: left_state.density,
                velocity: left_state.velocity,
                pressure: left_state.pressure,
            },
            DataPoint {
                coord: self.extrapolate(time),
                density: left_state.density,
                velocity: left_state.velocity,
                pressure: left_state.pressure,
            },
            DataPoint {
                coord: self.extrapolate(time) + shift,
                density: right_state.density,
                velocity: right_state.velocity,
                pressure: right_state.pressure,
            },
            DataPoint {
                coord: right_state.coord,
                density: right_state.density,
                velocity: right_state.velocity,
                pressure: right_state.pressure,
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
                        density: state.density,
                        velocity: state.velocity,
                        pressure: state.pressure,
                    },
                    // contact surface data
                    DataPoint {
                        coord: shock_coord + shift,
                        density: self.density(state),
                        velocity: contact.velocity,
                        pressure: self.pressure,
                    },
                ]
            }
            StateSide::Right => {
                vec![
                    // contact surface data first
                    DataPoint {
                        coord: shock_coord,
                        density: self.density(state),
                        velocity: contact.velocity,
                        pressure: self.pressure,
                    },
                    // right state data
                    DataPoint {
                        coord: shock_coord + shift,
                        density: state.density,
                        velocity: state.velocity,
                        pressure: state.pressure,
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

    /// Calculate adjusted density in a shocked state.
    /// Typically done after the iterative process is over.
    pub fn density(&self, state: &EulerState) -> f64 {
        state.gamma * self.pressure / (self.sound_speed(state) * self.sound_speed(state))
    }

    /// Calculate adjusted sound speed because of a shock.
    /// Typically done after the iterative process is over.
    pub fn sound_speed(&self, state: &EulerState) -> f64 {
        let top_term = state.gamma + 1.0 + (state.gamma - 1.0) * self.pressure / state.pressure;
        let bottom_term = state.gamma + 1.0 + (state.gamma - 1.0) * state.pressure / self.pressure;

        state.sound_speed() * (top_term / bottom_term).sqrt()
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

        // choose 12 points to plot -- arbitrary
        let interpolation_points = 20;
        // + head, + tail
        let total_points = 1 + interpolation_points + 1;

        let lower_vel = state.velocity.min(contact.velocity);
        let higher_vel = state.velocity.max(contact.velocity);

        let velocity_step = (higher_vel - lower_vel) / total_points as f64;
        let mut curr_vel = lower_vel;

        let mut raref_data: Vec<_> = Vec::with_capacity(total_points);

        for _ in 0..total_points+1 {

            raref_data.push(
                DataPoint {
                    coord: Self::extrapolate(state, curr_vel, Self::sound_speed(state, curr_vel), time),
                    density: Self::density(state, curr_vel),
                    velocity: curr_vel,
                    pressure: Self::pressure(state, Self::sound_speed(state, curr_vel)),
            });

            curr_vel += velocity_step;
            if curr_vel > higher_vel {
                curr_vel = higher_vel;
            }
        }

        let head_pos = Self::extrapolate(state, state.velocity, Self::sound_speed(state, state.velocity), time);
        let tail_pos = Self::extrapolate(state, contact.velocity, Self::sound_speed(state, contact.velocity), time);

        if head_pos > tail_pos {
            raref_data.reverse()
        }

        raref_data
    }

    /// Extrapolate a rarefaction wave given the initial state and a contact surface.
    pub fn extrapolate(state: &EulerState, velocity: f64, sound_speed: f64, time: f64) -> f64 {
         match state.side {
            StateSide::Left => {
                (velocity - sound_speed) * time
            }
            StateSide::Right => {
                (velocity + sound_speed) * time
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

    /// Rarefaction pressure for a given sound speed.
    fn pressure(state: &EulerState, sound_speed: f64) -> f64 {
        let gamma_expo = 2.0 * state.gamma / (state.gamma - 1.0);
        state.pressure * (sound_speed / state.sound_speed()).powf(gamma_expo)
    }

    /// Rarefaction pressure derivative relative to velocity for a given sound speed.
    fn pressure_derivative(state: &EulerState, sound_speed: f64, pressure: f64) -> f64 {
        match state.side {
            StateSide::Left => -state.gamma * pressure / sound_speed,
            StateSide::Right => state.gamma * pressure / sound_speed,
        }
    }
}
