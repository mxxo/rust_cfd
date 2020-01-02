// ----------------------------------------------------------------------------
//  Dependencies
// ----------------------------------------------------------------------------

// floating-point equality
extern crate approx;

use crate::fluxes::*;
use crate::limiters::VanAlbada;
use crate::pde::Boundary1d;
use crate::riemann::{waves::DataPoint, DomainBounds, EulerState, StateSide};

use std::ops::{Add, Mul, Sub};

// the solution vector
#[derive(Debug, Clone)]
pub struct EulerSolution1d {
    // cell values
    cells: Vec<EulerCell1d>,
    // the discretized domain
    domain: Vec<Boundary1d>,
    // the smallest cell width
    pub delta_x: f64,
    pub cfl: f64,
}

impl EulerSolution1d {
    pub fn init(
        left: PrimitiveState,
        right: PrimitiveState,
        bounds: DomainBounds,
        num_cells: usize,
        cfl: f64,
    ) -> Self {
        // initialize the domain
        let domain = domain_from_bounds(bounds, num_cells);

        // reserve space for the cell vector
        let mut cells: Vec<EulerCell1d> = Vec::with_capacity(num_cells);

        // initialize solution vector
        for cell_idx in 0..num_cells {
            let left_coord = domain[cell_idx].coord;
            let right_coord = domain[cell_idx + 1].coord;

            let center = (right_coord + left_coord) / 2.0;

            if center <= bounds.interface {
                cells.push(EulerCell1d::new(left, (left_coord, right_coord)));
            } else {
                cells.push(EulerCell1d::new(right, (left_coord, right_coord)));
            }
        }

        let delta_x = cells
            .iter()
            .map(|cell| cell.right - cell.left) // extract cell widths
            .fold(std::f64::INFINITY, |a, b| a.min(b)); // reduce to minimum value

        Self {
            cells,
            domain,
            delta_x,
            cfl,
        }
    }

    // max (absolute) wave speed is either
    //   |u - a| or |u + a|
    // where u can be negative. we can simplify this
    // by finding |u| + a since the sound speed is always positive
    pub fn max_stable_timestep(&self) -> f64 {
        self.cfl * self.delta_x / self.max_wave_speed()
    }

    pub fn max_wave_speed(&self) -> f64 {
        self.cells
            .iter()
            .map(|cell| cell.velocity().abs() + cell.sound_speed())
            .fold(0.0, |a, b| a.max(b)) // reduce to maximum
    }

    /// Time march this solution until a final time.
    pub fn first_order_time_march(
        mut self,
        flux_fn: impl FluxFunction,
        t_final: f64,
    ) -> Vec<EulerCell1d> {
        let mut t = 0.0;
        while t < t_final {
            // find
            // skim off previous values
            let old_soln = self.clone();
            let time_step = if t + old_soln.max_stable_timestep() > t_final {
                t_final - t
            } else {
                old_soln.max_stable_timestep()
            };

            // ignore boundaries for now :)), we don't really care about them
            // for our four reference problems
            for i in 1..old_soln.cells.len() - 1 {
                let left_flux =
                    flux_fn.calculate_flux(old_soln.cells[i - 1], old_soln.cells[i], time_step);
                let right_flux =
                    flux_fn.calculate_flux(old_soln.cells[i], old_soln.cells[i + 1], time_step);

                self.cells[i].advance(time_step, left_flux, right_flux);
            }

            t += time_step;
        }

        self.cells
    }

    /// Second order time-march with predictor-corrector and linear reconstruction.
    pub fn second_order_time_march(
        mut self,
        flux_fn: impl FluxFunction,
        t_final: f64,
    ) -> Vec<EulerCell1d> {
        let mut t = 0.0;
        while t < t_final {
            // skim off previous values
            let old_soln = self.clone();

            let time_step = if t + old_soln.max_stable_timestep() > t_final {
                t_final - t
            } else {
                old_soln.max_stable_timestep()
            };

            // -- prediction

            // find old solution slope limits
            let old_limits = VanAlbada::soln_limiters(&old_soln.cells);

            // fill a vector with all the guess values
            let guesses = {
                let mut guesses = old_soln.cells.clone();

                // don't update first and last cells for simplicity
                for i in 1..guesses.len() - 1 {
                    let left_flux = flux_fn.calculate_flux(
                        old_soln.cells[i - 1].reconstruct(old_limits[i - 1]),
                        old_soln.cells[i].reconstruct(-1.0 * old_limits[i]),
                        time_step,
                    );
                    let right_flux = flux_fn.calculate_flux(
                        old_soln.cells[i].reconstruct(old_limits[i]),
                        old_soln.cells[i + 1].reconstruct(-1.0 * old_limits[i + 1]),
                        time_step,
                    );

                    guesses[i] = old_soln.cells[i].cell_guess(time_step, left_flux, right_flux);
                }

                guesses
            };

            // find predicted solution limits
            let guess_limits = VanAlbada::soln_limiters(&guesses);

            // dbg!(guesses);

            // -- correction

            // use guesses to adjust cell updates
            for i in 1..old_soln.cells.len() - 1 {
                // let left_flux =
                //     flux_fn.calculate_flux(old_soln.cells[i - 1], old_soln.cells[i], time_step);
                // let right_flux =
                //     flux_fn.calculate_flux(old_soln.cells[i], old_soln.cells[i + 1], time_step);

                // old values
                let left_flux = flux_fn.calculate_flux(
                    old_soln.cells[i - 1].reconstruct(old_limits[i - 1]),
                    old_soln.cells[i].reconstruct(-1.0 * old_limits[i]),
                    time_step,
                );
                let right_flux = flux_fn.calculate_flux(
                    old_soln.cells[i].reconstruct(old_limits[i]),
                    old_soln.cells[i + 1].reconstruct(-1.0 * old_limits[i + 1]),
                    time_step,
                );

                // guess values
                let left_flux_guess = flux_fn.calculate_flux(
                    guesses[i - 1].reconstruct(guess_limits[i - 1]),
                    guesses[i].reconstruct(-1.0 * guess_limits[i]),
                    time_step,
                );

                let right_flux_guess = flux_fn.calculate_flux(
                    guesses[i].reconstruct(guess_limits[i]),
                    guesses[i + 1].reconstruct(-1.0 * guess_limits[i + 1]),
                    time_step,
                );

                // update cell with guess-averaged values
                // -- combine left, right fluxes
                let left_avg = 0.5 * (left_flux + left_flux_guess);
                let right_avg = 0.5 * (right_flux + right_flux_guess);

                self.cells[i].advance(time_step, left_avg, right_avg);
            }

            t += time_step;
        }

        self.cells
    }
}

fn domain_from_bounds(bounds: DomainBounds, num_cells: usize) -> Vec<Boundary1d> {
    use crate::pde::make_domain;
    make_domain(bounds.left, bounds.right, num_cells)
}

/// A primitive set of variables for the Euler equations.
#[derive(Debug, Clone, Copy)]
pub struct PrimitiveState {
    pub density: f64,
    pub velocity: f64,
    pub pressure: f64,
    pub gamma: f64,
}

impl PrimitiveState {
    pub fn from_data_point(val: DataPoint, gamma: f64) -> Self {
        Self {
            density: val.density,
            velocity: val.velocity,
            pressure: val.pressure,
            gamma,
        }
    }

    pub fn sound_speed(&self) -> f64 {
        (self.gamma * self.pressure / self.density).sqrt()
    }

    pub fn to_flux(self) -> EulerFlux {
        EulerFlux::new(self)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct PrimitiveResult {
    pub state: PrimitiveState,
    pub coord: f64,
    pub width: f64,
}

/// A discrete Euler 1D solution cell.
#[derive(Debug, Clone, Copy)]
pub struct EulerCell1d {
    pub density: f64,
    pub momentum: f64,
    pub energy: f64,
    pub gamma: f64,
    pub left: f64,
    pub right: f64,
}

impl EulerCell1d {
    /// Advance this cell forward in time.
    pub fn advance(&mut self, timestep: f64, left_flux: EulerFlux, right_flux: EulerFlux) {
        // update this cell
        *self = self.cell_guess(timestep, left_flux, right_flux) // + scaling_factor * net_flux;
    }

    /// Get a copy of an updated cell.
    pub fn cell_guess(
        self,
        timestep: f64,
        left_flux: EulerFlux,
        right_flux: EulerFlux,
    ) -> EulerCell1d {
        let net_flux = left_flux - right_flux;
        let scaling_factor = timestep / self.width();

        self + scaling_factor * net_flux
    }

    /// Make a new Euler state from primitive variables.
    pub fn new(state: PrimitiveState, bounds: (f64, f64)) -> Self {
        let (left, right) = bounds;
        Self {
            density: state.density,
            momentum: state.density * state.velocity,
            energy: Self::energy(state),
            gamma: state.gamma,
            left,
            right,
        }
    }

    // -- state equations
    pub fn to_primitive_result(self) -> PrimitiveResult {
        PrimitiveResult {
            width: self.width(),
            coord: self.center(),
            state: self.to_primitive(),
        }
    }

    pub fn to_primitive(self) -> PrimitiveState {
        PrimitiveState {
            density: self.density,
            velocity: self.velocity(),
            pressure: self.pressure(),
            gamma: self.gamma,
        }
    }

    pub fn to_flux(self) -> EulerFlux {
        self.to_primitive().to_flux()
    }

    pub fn to_euler_state(self, side: StateSide) -> EulerState {
        EulerState {
            density: self.density,
            velocity: self.velocity(),
            pressure: self.pressure(),
            gamma: self.gamma,
            side,
        }
    }

    pub fn velocity(&self) -> f64 {
        self.momentum / self.density
    }

    pub fn pressure(&self) -> f64 {
        (self.gamma - 1.0)
            * (self.energy - (self.density * self.velocity() * self.velocity() / 2.0))
    }

    pub fn sound_speed(&self) -> f64 {
        (self.gamma * self.pressure() / self.density).sqrt()
    }

    /// Energy expression for Euler.
    fn energy(state: PrimitiveState) -> f64 {
        state.pressure / (state.gamma - 1.0)
            + (state.density * state.velocity * state.velocity / 2.0)
    }

    /// Cell center.
    pub fn center(&self) -> f64 {
        (self.right + self.left) / 2.0
    }

    /// Cell width.
    pub fn width(&self) -> f64 {
        self.right - self.left
    }
}

impl Add<EulerFlux> for EulerCell1d {
    type Output = Self;

    fn add(self, rhs: EulerFlux) -> Self::Output {
        Self {
            density: self.density + rhs.density_flux,
            momentum: self.momentum + rhs.momentum_flux,
            energy: self.energy + rhs.energy_flux,
            gamma: self.gamma,
            left: self.left,
            right: self.right,
        }
    }
}

/// An Euler flux set.
#[derive(Debug, Clone, Copy)]
pub struct EulerFlux {
    pub density_flux: f64,
    pub momentum_flux: f64,
    pub energy_flux: f64,
}

impl EulerFlux {
    /// Get a new flux vector based on the primitive state.
    pub fn new(state: PrimitiveState) -> Self {
        Self {
            density_flux: state.density * state.velocity,
            momentum_flux: Self::momentum_flux(state),
            energy_flux: Self::energy_flux(state),
        }
    }

    /// Calculate momentum flux based on primitive variables.
    fn momentum_flux(state: PrimitiveState) -> f64 {
        state.density * state.velocity * state.velocity + state.pressure
    }

    /// Calculate energy flux based on primitive variables.
    fn energy_flux(state: PrimitiveState) -> f64 {
        state.velocity
            * (state.gamma * state.pressure / (state.gamma - 1.0)
                + (state.density * state.velocity * state.velocity / 2.0))
    }
}

// net flux
impl Sub for EulerFlux {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            density_flux: self.density_flux - other.density_flux,
            momentum_flux: self.momentum_flux - other.momentum_flux,
            energy_flux: self.energy_flux - other.energy_flux,
        }
    }
}

impl Add<EulerFlux> for EulerFlux {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self {
            density_flux: self.density_flux + rhs.density_flux,
            momentum_flux: self.momentum_flux + rhs.momentum_flux,
            energy_flux: self.energy_flux + rhs.energy_flux,
        }
    }
}

// scalar multiplication
impl Mul<EulerFlux> for f64 {
    type Output = EulerFlux;

    fn mul(self, rhs: EulerFlux) -> Self::Output {
        EulerFlux {
            density_flux: self * rhs.density_flux,
            momentum_flux: self * rhs.momentum_flux,
            energy_flux: self * rhs.energy_flux,
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use approx::*;

    #[test]
    fn add_multiply_flux() {
        let left_state = PrimitiveState {
            density: 2.281,
            velocity: 164.83,
            pressure: 201.17e3,
            gamma: 1.4,
        };

        let double_state = left_state.to_flux() + left_state.clone().to_flux();
        assert_relative_eq!(
            (2.0 * left_state.to_flux()).density_flux,
            double_state.density_flux,
            epsilon = 0.001
        );
    }

    // #[test]
    // fn subtract_euler_state() {}
}
