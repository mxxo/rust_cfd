// ----------------------------------------------------------------------------
//  Dependencies
// ----------------------------------------------------------------------------

// floating-point equality
extern crate approx;

use crate::riemann::DomainBounds;
use std::ops::{Mul, Sub};

// the solution vector
#[derive(Debug)]
pub struct EulerSolution1d {
    // cell values
    cells: Vec<EulerCell1d>,
    // the discretized domain
    domain: Vec<Boundary1d>,
    // the smallest cell width
    pub delta_x: f64,
}


impl EulerSolution1d {

    pub fn init(left: PrimitiveState, right: PrimitiveState, bounds: DomainBounds, num_cells: u32) -> Self {
        unimplemented!();
        // Self {
        //
        // }
    }
}

// ----------------------------------------------------------------------------
// Flux functions
// ----------------------------------------------------------------------------

/// Possible flux functions for Euler.
#[derive(Debug, Clone, Copy)]
pub enum FluxType {
    Exact,
    HLLE,
    Roe,
    RoeEntropy,
}

/// A primitive set of variables for the Euler equations.
#[derive(Debug, Clone, Copy)]
pub struct PrimitiveState {
    pub density: f64,
    pub velocity: f64,
    pub pressure: f64,
    pub gamma: f64,
}

/// A discrete Euler 1D solution cell.
#[derive(Debug, Clone, Copy)]
pub struct EulerCell1d {
    pub density: f64,
    pub momentum: f64,
    pub energy: f64,
    pub left: f64,
    pub right: f64,
}

impl EulerCell1d {

    /// Advance this cell forward in time.
    pub fn advance(&mut self, timestep: f64, left_flux: EulerFlux, right_flux: EulerFlux) {
        // flux entering - flux leaving
        let net_flux = left_flux - right_flux;
        let scaling_factor = timestep / (self.right - self.left);
        let (new_dens, new_mom, new_energ) = net_flux * scaling_factor;
        self.density = new_dens;
        self.momentum = new_mom;
        self.energy = new_energ;
    }


    /// Make a new Euler state from primitive variables.
    pub fn new(state: PrimitiveState, bounds: (f64, f64)) -> Self {
        let (left, right) = bounds;
        Self {
            density: state.density,
            momentum: state.density * state.velocity,
            energy: Self::energy(state),
            left,
            right,
        }
    }

    /// Energy expression for Euler.
    fn energy(state: PrimitiveState) -> f64 {
        state.pressure / (state.gamma - 1.0) + (state.density * state.velocity * state.velocity / 2.0)
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
            momentum_flux: Self::get_momentum_flux(state),
            energy_flux: Self::get_energy_flux(state),
        }
    }

    /// Calculate momentum flux based on primitive variables.
    fn get_momentum_flux(state: PrimitiveState) -> f64 {
        state.density * state.velocity * state.velocity + state.pressure
    }

    /// Calculate energy flux based on primitive variables.
    fn get_energy_flux(state: PrimitiveState) -> f64 {
        state.velocity * (state.gamma * state.pressure / (state.gamma - 1.0) + (state.density * state.velocity * state.velocity / 2.0))
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

impl Mul<f64> for EulerFlux {
    type Output = (f64, f64, f64);

    fn mul(self, rhs: f64) -> Self::Output {
        (
            self.density_flux * rhs,
            self.momentum_flux * rhs,
            self.energy_flux * rhs,
        )
    }
}

// impl Cell1d {
//     // step this cell forward
//     pub fn advance(&mut self, timestep: f64, net_flux: f64) {
//         // some weirdness here since we're not using solution's uniform delta x?
//         self.mass += timestep / (self.right - self.left) * net_flux;
//     }
// }

// a boundary, aka cell interface
#[derive(Debug)]
struct Boundary1d {
    pub coord: f64,
    pub left_cell: usize,
    pub right_cell: usize,
}



//
// // flux based on the left cell only
// fn left_flux(boundary: &Boundary1d, past_values: &[f64], wavespeed: f64) -> f64 {
//     wavespeed * past_values[boundary.left_cell]
// }
//
// // flux based on the right cell only
// fn right_flux(boundary: &Boundary1d, past_values: &[f64], wavespeed: f64) -> f64 {
//     wavespeed * past_values[boundary.right_cell]
// }
//
// // F(average)
// fn average_flux(boundary: &Boundary1d, past_values: &[f64], wavespeed: f64) -> f64 {
//     wavespeed * 0.5 * (past_values[boundary.left_cell] + past_values[boundary.right_cell])
// }
//
// // F(Riemann solution)
// fn simple_riemann_flux(boundary: &Boundary1d, past_values: &[f64], wavespeed: f64) -> f64 {
//     // solution moving to the right, use information from the left cell
//     if wavespeed > 0.0 {
//         left_flux(boundary, past_values, wavespeed)
//     }
//     // solution moving to the left, use information from the right cell
//     else {
//         right_flux(boundary, past_values, wavespeed)
//     }
// }
//
// // ----------------------------------------------------------------------------
// // Method implementations
// // ----------------------------------------------------------------------------
//
// impl Solution1d {
//     // constructor
//     pub fn new(left: f64, right: f64, num_cells: usize) -> Self {
//         // initialize the domain
//         let domain = make_domain(left, right, num_cells);
//
//         // reserve space for the cell vector
//         let mut cells = Vec::with_capacity(num_cells);
//
//         // initialize solution vector
//         for cell_idx in 0..num_cells {
//             cells.push(Cell1d {
//                 value: 0.0,
//                 left: domain[cell_idx].coord,
//                 right: domain[cell_idx + 1].coord,
//             });
//         }
//
//         // find the minimum cell width
//         let delta_x = cells
//             .iter()
//             .map(|cell| cell.right - cell.left) // extract cell widths
//             .fold(std::f64::INFINITY, |a, b| a.min(b)); // reduce to minimum value
//
//         Solution1d {
//             cells,
//             domain,
//             delta_x,
//         }
//     }
//
//     // initialize the solution vector with an expression
//     pub fn init<F>(&mut self, ic_fn: F)
//     where
//         F: Fn(f64) -> f64,
//     {
//         // hardcode averaging scheme for simplicity
//         let quad = Midpoint::init(1000);
//         // average initial condition in each cell
//         for mut cell in self.cells.iter_mut() {
//             cell.value = quad.integrate(cell.left, cell.right, &ic_fn) / (cell.right - cell.left);
//         }
//     }
//
//     // fill the solution cell values from a slice
//     pub fn fill(&mut self, values: &[f64]) {
//         assert_eq!(values.len(), self.cells.len());
//
//         for (i, mut cell) in self.cells.iter_mut().enumerate() {
//             cell.value = values[i];
//         }
//     }
//
//     // update the solution
//     pub fn update(&mut self, flux_type: &FluxType, timestep: f64, wavespeed: f64) {
//         // use enum to dispatch the different flux functions
//         let flux_fn = match flux_type {
//             FluxType::LeftCell => left_flux,
//             FluxType::RightCell => right_flux,
//             FluxType::CellAverage => average_flux,
//             FluxType::Riemann => simple_riemann_flux,
//         };
//
//         // skim off previous timestep's values
//         // -- allocation seems to be factored out this method by optimizer
//         let past_values: Vec<f64> = self.cells.iter().map(|cell| cell.value).collect();
//
//         // update cells using the flux function
//         for (index, cell) in self.cells.iter_mut().enumerate() {
//             let left_boundary = &self.domain[index];
//             let right_boundary = &self.domain[index + 1];
//
//             let net_flux = flux_fn(left_boundary, &past_values, wavespeed)
//                 - flux_fn(right_boundary, &past_values, wavespeed);
//
//             cell.advance(timestep, net_flux);
//         }
//     }
//
//     // export the solution in a simple format, akin to csv
//     pub fn export(&self) -> (Vec<f64>, Vec<f64>) {
//         let coords = self
//             .cells
//             .iter()
//             .map(|cell| 0.5 * (cell.left + cell.right))
//             .collect();
//         let values = self.cells.iter().map(|cell| cell.value).collect();
//         (coords, values)
//     }
// }
//
// // helper function for making a domain
// //
// // create `num_cell` equally sized pieces between left and right
// fn make_domain(left: f64, right: f64, num_cells: usize) -> Vec<Boundary1d> {
//     // one more boundary than number of cells
//     let boundary_coords = linspace::<f64>(left, right, num_cells + 1);
//     let mut boundaries = Vec::new();
//
//     for (idx, coord) in boundary_coords.enumerate() {
//         // num_cells is for wraparound indexing of usize values
//         let prev_idx = if idx == 0 { num_cells - 1 } else { idx - 1 };
//         let next_idx = idx % num_cells;
//         boundaries.push(Boundary1d {
//             coord,
//             left_cell: prev_idx,
//             right_cell: next_idx,
//         });
//     }
//
//     boundaries
// }
//
// #[cfg(test)]
// mod tests {
//     use super::*;
//
//     #[test]
//     pub fn domain_size() {
//         let l = 100;
//         let domain = make_domain(-10., 10., l);
//         assert!(domain.len() == l + 1);
//     }
// }
