//! Assignment 5 - 2D Euler solver
// Max Orok January 2020

// ----------------------------------------------------------------------------
//  Dependencies
// ----------------------------------------------------------------------------
extern crate approx; // floating-point equality
extern crate num_traits; // safe unsigned-float conversions

// ----------------------------------------------------------------------------
// Library imports
// ----------------------------------------------------------------------------
use std::convert::From;
use std::ops::{Add, Mul, Sub};

use crate::euler1d::*;
use crate::fluxes::*;
use crate::limiters::VanAlbada;
use crate::pde::Boundary1d;
use crate::riemann::{waves::DataPoint, DomainBounds, EulerState, StateSide};

#[derive(Debug, Clone, Copy)]
pub struct Point2d {
    pub x: f64,
    pub y: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct EulerPrimitive2d {
    pub density: f64,
    pub x_vel: f64,
    pub y_vel: f64,
    pub pressure: f64,
    pub gamma: f64,
}

impl EulerPrimitive2d {
    pub fn energy(self) -> f64 {
        self.pressure / (self.gamma - 1.0)
            + (0.5 * self.density * (self.x_vel * self.x_vel + self.y_vel * self.y_vel))
    }
}

impl From<EulerPrimitive2d> for StateVec2d {
    fn from(prim: EulerPrimitive2d) -> Self {
        Self {
            density: prim.density,
            x_momentum: prim.density * prim.x_vel,
            y_momentum: prim.density * prim.y_vel,
            energy: prim.energy(),
            gamma: prim.gamma,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct StateVec2d {
    pub density: f64,
    pub x_momentum: f64,
    pub y_momentum: f64,
    pub energy: f64,
    pub gamma: f64,
}

/// An axis-aligned rectangular cell.
#[derive(Debug, Clone, Copy)]
pub struct EulerCell2d {
    pub state: StateVec2d,
    pub min: Point2d,
    pub max: Point2d,
}

impl EulerCell2d {
    pub fn centroid(&self) -> Point2d {
        Point2d {
            x: (self.min.x + self.max.x) / 2.0,
            y: (self.min.y + self.max.y) / 2.0,
        }
    }

    pub fn new(min: Point2d, max: Point2d) -> Self {
        Self {
            state: StateVec2d {
                density: 0.0,
                x_momentum: 0.0,
                y_momentum: 0.0,
                energy: 0.0,
                gamma: 0.0,
            },
            min,
            max,
        }
    }

    pub fn fill_state<T>(&mut self, new_state: T)
    where
        T: Into<StateVec2d>,
    {
        let new_state: StateVec2d = new_state.into();
        self.state = new_state;
    }
}

impl EulerCell2d {
    pub fn to_x_flux(self) -> EulerFlux2d {
        unimplemented!();
    }

    pub fn to_y_flux(self) -> EulerFlux2d {
        unimplemented!();
    }
}

/// A 2d Euler equations flux vector.
#[derive(Debug, Clone, Copy)]
pub struct EulerFlux2d {
    pub density: f64,
    pub x_momentum: f64,
    pub y_momentum: f64,
    pub energy: f64,
}

/// A grid of rectangular cells.
/// Could be extended to other shapes by replacing cells with
/// `Vec<Box<dyn EulerCell2d>>` for some trait `EulerCell2D`.
#[derive(Debug, Clone)]
pub struct EulerSolution2d {
    /// Cell states
    cells: Vec<EulerCell2d>,
    /// Number of x cells.
    pub num_x: usize,
    /// Number of y cells.
    pub num_y: usize,
    // the smallest cell boundaries
    pub delta_x: f64,
    pub delta_y: f64,
}

impl EulerSolution2d {
    /// Make an empty solution grid of `num_x` by `num_y` cells on the extents
    /// `min_point` to `max_point`.
    pub fn new(num_x: usize, num_y: usize, min_point: Point2d, max_point: Point2d) -> Self {
        use num_traits::cast::FromPrimitive;

        let mut cells: Vec<EulerCell2d> = Vec::with_capacity(num_x * num_y);

        let x_step = (max_point.x - min_point.x) / f64::from_usize(num_x).unwrap();
        let y_step = (max_point.y - min_point.y) / f64::from_usize(num_y).unwrap();

        for j in 0..num_y {
            for i in 0..num_x {
                cells.push(EulerCell2d::new(
                    Point2d {
                        x: x_step * f64::from_usize(i).unwrap() + min_point.x,
                        y: y_step * f64::from_usize(j).unwrap() + min_point.y,
                    },
                    Point2d {
                        x: x_step * f64::from_usize(i + 1).unwrap() + min_point.x,
                        y: y_step * f64::from_usize(j + 1).unwrap() + min_point.y,
                    },
                ));
            }
        }

        // find minimum cell boundaries

        let delta_x = cells
            .iter()
            .map(|cell| cell.max.x - cell.min.x) // extract cell widths
            .fold(std::f64::INFINITY, |a, b| a.min(b)); // reduce to minimum value

        let delta_y = cells
            .iter()
            .map(|cell| cell.max.y - cell.min.y) // extract cell widths
            .fold(std::f64::INFINITY, |a, b| a.min(b)); // reduce to minimum value

        Self {
            cells,
            num_x,
            num_y,
            delta_x,
            delta_y,
        }
    }

    /// Make an empty square solution grid of `side_cells` by `side_cells`.
    pub fn square(cells_width: usize, min_point: Point2d, max_point: Point2d) -> Self {
        Self::new(cells_width, cells_width, min_point, max_point)
    }

    /// Populate the grid using an expression based on each cell's `x-y` position.
    pub fn init<F, T>(&mut self, pos_expr: F)
    where
        F: Fn(Point2d) -> T,
        T: Into<StateVec2d>,
    {
        for cell in &mut self.cells {
            cell.fill_state(pos_expr(cell.centroid().into()));
        }
    }

    /// Get the cell at `[x_idx][y_idx]`. Uses zero-indexing.
    pub fn at(&self, x_idx: usize, y_idx: usize) -> EulerCell2d {
        // built-in bounds-checking
        self.cells[self.index(x_idx, y_idx)]
    }

    #[inline(always)]
    fn index(&self, x_idx: usize, y_idx: usize) -> usize {
        x_idx + self.num_x * y_idx
    }

    // solver functions
    /// Time march this solution until a final time.
    pub fn first_order_time_march(
        mut self,
        cfl: f64,
        flux_fn: impl FluxFunction,
        t_final: f64,
    ) -> Vec<EulerCell2d> {
        self.cells

        // let mut t = 0.0;
        // while t < t_final {
        //     // find
        //     // skim off previous values
        //     let old_soln = self.clone();
        //     let time_step = if t + old_soln.max_stable_timestep() > t_final {
        //         t_final - t
        //     } else {
        //         old_soln.max_stable_timestep()
        //     };

        //     // ignore boundaries for now :)), we don't really care about them
        //     // for our four reference problems
        //     for i in 1..old_soln.cells.len() - 1 {
        //         let left_flux =
        //             flux_fn.calculate_flux(old_soln.cells[i - 1], old_soln.cells[i], time_step);
        //         let right_flux =
        //             flux_fn.calculate_flux(old_soln.cells[i], old_soln.cells[i + 1], time_step);

        //         self.cells[i].advance(time_step, left_flux, right_flux);
        //     }

        //     t += time_step;
        // }

        // self.cells
    }

    // data viz for cells
    pub fn nodes(&self) -> Vec<Point2d> {
        let mut result = Vec::with_capacity(self.cells.len() + self.num_x + 1);
        dbg!(result.len());

        // scrape off the bottom nodes of each row
        for (idx, cell) in self.cells.iter().enumerate() {
            // first column of cells need their min cell values
            if idx % self.num_x == 0 {
                result.push(Point2d {
                    x: cell.min.x,
                    y: cell.min.y,
                });
            }
            // always push this cell
            result.push(Point2d {
                x: cell.max.x,
                y: cell.min.y,
            });
        }

        // take top points of top row of cells
        for (idx, cell) in self.cells[self.cells.len() - self.num_x..]
            .iter()
            .enumerate()
        {
            if idx % self.num_x == 0 {
                result.push(Point2d {
                    x: cell.min.x,
                    y: cell.max.y,
                });
            }
            // always push this cell
            result.push(Point2d {
                x: cell.max.x,
                y: cell.max.y,
            });
        }

        result
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use approx::*;

    #[test]
    fn indexing() {
        let soln =
            EulerSolution2d::square(2, Point2d { x: -0.5, y: -0.5 }, Point2d { x: 0.5, y: 0.5 });

        let cell = soln.at(0, 0);

        assert_relative_eq!(cell.centroid().x, -0.25);
        assert_relative_eq!(cell.centroid().y, -0.25);

        let cell = soln.at(1, 0);
        assert_relative_eq!(cell.centroid().x, 0.25);
        assert_relative_eq!(cell.centroid().y, -0.25);

        let cell = soln.at(0, 1);
        assert_relative_eq!(cell.centroid().x, -0.25);
        assert_relative_eq!(cell.centroid().y, 0.25);

        let cell = soln.at(1, 1);
        assert_relative_eq!(cell.centroid().x, 0.25);
        assert_relative_eq!(cell.centroid().y, 0.25);
    }
}

// impl EulerSolution1d {
//     pub fn init(
//         left: PrimitiveState,
//         right: PrimitiveState,
//         bounds: DomainBounds,
//         num_cells: usize,
//         cfl: f64,
//     ) -> Self {
//         // initialize the domain
//         let domain = domain_from_bounds(bounds, num_cells);
//
//         // reserve space for the cell vector
//         let mut cells: Vec<EulerCell1d> = Vec::with_capacity(num_cells);
//
//         // initialize solution vector
//         for cell_idx in 0..num_cells {
//             let left_coord = domain[cell_idx].coord;
//             let right_coord = domain[cell_idx + 1].coord;
//
//             let center = (right_coord + left_coord) / 2.0;
//
//             if center <= bounds.interface {
//                 cells.push(EulerCell1d::new(left, (left_coord, right_coord)));
//             } else {
//                 cells.push(EulerCell1d::new(right, (left_coord, right_coord)));
//             }
//         }
//
//         let delta_x = cells
//             .iter()
//             .map(|cell| cell.right - cell.left) // extract cell widths
//             .fold(std::f64::INFINITY, |a, b| a.min(b)); // reduce to minimum value
//
//         Self {
//             cells,
//             domain,
//             delta_x,
//             cfl,
//         }
//     }
//
//     // max (absolute) wave speed is either
//     //   |u - a| or |u + a|
//     // where u can be negative. we can simplify this
//     // by finding |u| + a since the sound speed is always positive
//     pub fn max_stable_timestep(&self) -> f64 {
//         self.cfl * self.delta_x / self.max_wave_speed()
//     }
//
//     pub fn max_wave_speed(&self) -> f64 {
//         self.cells
//             .iter()
//             .map(|cell| cell.velocity().abs() + cell.sound_speed())
//             .fold(0.0, |a, b| a.max(b)) // reduce to maximum
//     }
//
//     /// Time march this solution until a final time.
//     pub fn first_order_time_march(
//         mut self,
//         flux_fn: impl FluxFunction,
//         t_final: f64,
//     ) -> Vec<EulerCell1d> {
//         let mut t = 0.0;
//         while t < t_final {
//             // find
//             // skim off previous values
//             let old_soln = self.clone();
//             let time_step = if t + old_soln.max_stable_timestep() > t_final {
//                 t_final - t
//             } else {
//                 old_soln.max_stable_timestep()
//             };
//
//             // ignore boundaries for now :)), we don't really care about them
//             // for our four reference problems
//             for i in 1..old_soln.cells.len() - 1 {
//                 let left_flux =
//                     flux_fn.calculate_flux(old_soln.cells[i - 1], old_soln.cells[i], time_step);
//                 let right_flux =
//                     flux_fn.calculate_flux(old_soln.cells[i], old_soln.cells[i + 1], time_step);
//
//                 self.cells[i].advance(time_step, left_flux, right_flux);
//             }
//
//             t += time_step;
//         }
//
//         self.cells
//     }
//
//     /// Second order time-march with predictor-corrector and linear reconstruction.
//     pub fn second_order_time_march(
//         mut self,
//         flux_fn: impl FluxFunction,
//         t_final: f64,
//     ) -> Vec<EulerCell1d> {
//         let mut t = 0.0;
//         while t < t_final {
//             // skim off previous values
//             let old_soln = self.clone();
//
//             let time_step = if t + old_soln.max_stable_timestep() > t_final {
//                 t_final - t
//             } else {
//                 old_soln.max_stable_timestep()
//             };
//
//             // -- prediction
//
//             // find old solution slope limits
//             let old_limits = VanAlbada::soln_limiters(&old_soln.cells);
//
//             // fill a vector with all the guess values
//             let guesses = {
//                 let mut guesses = old_soln.cells.clone();
//
//                 // don't update first and last cells for simplicity
//                 for i in 1..guesses.len() - 1 {
//                     let left_flux = flux_fn.calculate_flux(
//                         old_soln.cells[i - 1].reconstruct(old_limits[i - 1]),
//                         old_soln.cells[i].reconstruct(-1.0 * old_limits[i]),
//                         time_step,
//                     );
//                     let right_flux = flux_fn.calculate_flux(
//                         old_soln.cells[i].reconstruct(old_limits[i]),
//                         old_soln.cells[i + 1].reconstruct(-1.0 * old_limits[i + 1]),
//                         time_step,
//                     );
//
//                     guesses[i] = old_soln.cells[i].cell_guess(time_step, left_flux, right_flux);
//                 }
//
//                 guesses
//             };
//
//             // find predicted solution limits
//             let guess_limits = VanAlbada::soln_limiters(&guesses);
//
//             // dbg!(guesses);
//
//             // -- correction
//
//             // use guesses to adjust cell updates
//             for i in 1..old_soln.cells.len() - 1 {
//                 // let left_flux =
//                 //     flux_fn.calculate_flux(old_soln.cells[i - 1], old_soln.cells[i], time_step);
//                 // let right_flux =
//                 //     flux_fn.calculate_flux(old_soln.cells[i], old_soln.cells[i + 1], time_step);
//
//                 // old values
//                 let left_flux = flux_fn.calculate_flux(
//                     old_soln.cells[i - 1].reconstruct(old_limits[i - 1]),
//                     old_soln.cells[i].reconstruct(-1.0 * old_limits[i]),
//                     time_step,
//                 );
//                 let right_flux = flux_fn.calculate_flux(
//                     old_soln.cells[i].reconstruct(old_limits[i]),
//                     old_soln.cells[i + 1].reconstruct(-1.0 * old_limits[i + 1]),
//                     time_step,
//                 );
//
//                 // guess values
//                 let left_flux_guess = flux_fn.calculate_flux(
//                     guesses[i - 1].reconstruct(guess_limits[i - 1]),
//                     guesses[i].reconstruct(-1.0 * guess_limits[i]),
//                     time_step,
//                 );
//
//                 let right_flux_guess = flux_fn.calculate_flux(
//                     guesses[i].reconstruct(guess_limits[i]),
//                     guesses[i + 1].reconstruct(-1.0 * guess_limits[i + 1]),
//                     time_step,
//                 );
//
//                 // update cell with guess-averaged values
//                 // -- combine left, right fluxes
//                 let left_avg = 0.5 * (left_flux + left_flux_guess);
//                 let right_avg = 0.5 * (right_flux + right_flux_guess);
//
//                 self.cells[i].advance(time_step, left_avg, right_avg);
//             }
//
//             t += time_step;
//         }
//
//         self.cells
//     }
// }
//
// fn domain_from_bounds(bounds: DomainBounds, num_cells: usize) -> Vec<Boundary1d> {
//     use crate::pde::make_domain;
//     make_domain(bounds.left, bounds.right, num_cells)
// }
//
// /// A primitive set of variables for the Euler equations.
// #[derive(Debug, Clone, Copy)]
// pub struct PrimitiveState {
//     pub density: f64,
//     pub velocity: f64,
//     pub pressure: f64,
//     pub gamma: f64,
// }
//
// impl PrimitiveState {
//     pub fn from_data_point(val: DataPoint, gamma: f64) -> Self {
//         Self {
//             density: val.density,
//             velocity: val.velocity,
//             pressure: val.pressure,
//             gamma,
//         }
//     }
//
//     pub fn sound_speed(&self) -> f64 {
//         (self.gamma * self.pressure / self.density).sqrt()
//     }
//
//     pub fn to_flux(self) -> EulerFlux {
//         EulerFlux::new(self)
//     }
// }
//
// #[derive(Debug, Clone, Copy)]
// pub struct PrimitiveResult {
//     pub state: PrimitiveState,
//     pub coord: f64,
//     pub width: f64,
// }
//
// /// A discrete Euler 2D solution cell.
// #[derive(Debug, Clone, Copy)]
// pub struct EulerCell2d {
//     pub density: f64,
//     pub momentum: f64,
//     pub energy: f64,
//     pub gamma: f64,
//     pub left: f64,
//     pub right: f64,
// }
//
// impl EulerCell1d {
//     /// Advance this cell forward in time.
//     pub fn advance(&mut self, timestep: f64, left_flux: EulerFlux, right_flux: EulerFlux) {
//         // update this cell
//         *self = self.cell_guess(timestep, left_flux, right_flux) // + scaling_factor * net_flux;
//     }
//
//     /// Get a copy of an updated cell.
//     pub fn cell_guess(
//         self,
//         timestep: f64,
//         left_flux: EulerFlux,
//         right_flux: EulerFlux,
//     ) -> EulerCell1d {
//         let net_flux = left_flux - right_flux;
//         let scaling_factor = timestep / self.width();
//
//         self + scaling_factor * net_flux
//     }
//
//     /// Make a new Euler state from primitive variables.
//     pub fn new(state: PrimitiveState, bounds: (f64, f64)) -> Self {
//         let (left, right) = bounds;
//         Self {
//             density: state.density,
//             momentum: state.density * state.velocity,
//             energy: Self::energy(state),
//             gamma: state.gamma,
//             left,
//             right,
//         }
//     }
//
//     // -- state equations
//     pub fn to_primitive_result(self) -> PrimitiveResult {
//         PrimitiveResult {
//             width: self.width(),
//             coord: self.center(),
//             state: self.to_primitive(),
//         }
//     }
//
//     pub fn to_primitive(self) -> PrimitiveState {
//         PrimitiveState {
//             density: self.density,
//             velocity: self.velocity(),
//             pressure: self.pressure(),
//             gamma: self.gamma,
//         }
//     }
//
//     pub fn to_flux(self) -> EulerFlux {
//         self.to_primitive().to_flux()
//     }
//
//     pub fn to_euler_state(self, side: StateSide) -> EulerState {
//         EulerState {
//             density: self.density,
//             velocity: self.velocity(),
//             pressure: self.pressure(),
//             gamma: self.gamma,
//             side,
//         }
//     }
//
//     pub fn velocity(&self) -> f64 {
//         self.momentum / self.density
//     }
//
//     pub fn pressure(&self) -> f64 {
//         (self.gamma - 1.0)
//             * (self.energy - (self.density * self.velocity() * self.velocity() / 2.0))
//     }
//
//     pub fn sound_speed(&self) -> f64 {
//         (self.gamma * self.pressure() / self.density).sqrt()
//     }
//
//     /// Energy expression for Euler.
//     fn energy(state: PrimitiveState) -> f64 {
//         state.pressure / (state.gamma - 1.0)
//             + (state.density * state.velocity * state.velocity / 2.0)
//     }
//
//     /// Cell center.
//     pub fn center(&self) -> f64 {
//         (self.right + self.left) / 2.0
//     }
//
//     /// Cell width.
//     pub fn width(&self) -> f64 {
//         self.right - self.left
//     }
// }
//
// impl Add<EulerFlux> for EulerCell1d {
//     type Output = Self;
//
//     fn add(self, rhs: EulerFlux) -> Self::Output {
//         Self {
//             density: self.density + rhs.density_flux,
//             momentum: self.momentum + rhs.momentum_flux,
//             energy: self.energy + rhs.energy_flux,
//             gamma: self.gamma,
//             left: self.left,
//             right: self.right,
//         }
//     }
// }
//
// /// An Euler flux set.
// #[derive(Debug, Clone, Copy)]
// pub struct EulerFlux {
//     pub density_flux: f64,
//     pub momentum_flux: f64,
//     pub energy_flux: f64,
// }
//
// impl EulerFlux {
//     /// Get a new flux vector based on the primitive state.
//     pub fn new(state: PrimitiveState) -> Self {
//         Self {
//             density_flux: state.density * state.velocity,
//             momentum_flux: Self::momentum_flux(state),
//             energy_flux: Self::energy_flux(state),
//         }
//     }
//
//     /// Calculate momentum flux based on primitive variables.
//     fn momentum_flux(state: PrimitiveState) -> f64 {
//         state.density * state.velocity * state.velocity + state.pressure
//     }
//
//     /// Calculate energy flux based on primitive variables.
//     fn energy_flux(state: PrimitiveState) -> f64 {
//         state.velocity
//             * (state.gamma * state.pressure / (state.gamma - 1.0)
//                 + (state.density * state.velocity * state.velocity / 2.0))
//     }
// }
//
// // net flux
// impl Sub for EulerFlux {
//     type Output = Self;
//
//     fn sub(self, other: Self) -> Self {
//         Self {
//             density_flux: self.density_flux - other.density_flux,
//             momentum_flux: self.momentum_flux - other.momentum_flux,
//             energy_flux: self.energy_flux - other.energy_flux,
//         }
//     }
// }
//
// impl Add<EulerFlux> for EulerFlux {
//     type Output = Self;
//
//     fn add(self, rhs: Self) -> Self {
//         Self {
//             density_flux: self.density_flux + rhs.density_flux,
//             momentum_flux: self.momentum_flux + rhs.momentum_flux,
//             energy_flux: self.energy_flux + rhs.energy_flux,
//         }
//     }
// }
//
// // scalar multiplication
// impl Mul<EulerFlux> for f64 {
//     type Output = EulerFlux;
//
//     fn mul(self, rhs: EulerFlux) -> Self::Output {
//         EulerFlux {
//             density_flux: self * rhs.density_flux,
//             momentum_flux: self * rhs.momentum_flux,
//             energy_flux: self * rhs.energy_flux,
//         }
//     }
// }
//
// #[cfg(test)]
// mod tests {
//
//     use super::*;
//     use approx::*;
//
//     #[test]
//     fn add_multiply_flux() {
//         let left_state = PrimitiveState {
//             density: 2.281,
//             velocity: 164.83,
//             pressure: 201.17e3,
//             gamma: 1.4,
//         };
//
//         let double_state = left_state.to_flux() + left_state.clone().to_flux();
//         assert_relative_eq!(
//             (2.0 * left_state.to_flux()).density_flux,
//             double_state.density_flux,
//             epsilon = 0.001
//         );
//     }
//
//     // #[test]
//     // fn subtract_euler_state() {}
// }
