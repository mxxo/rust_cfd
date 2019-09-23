// Finite-Volume methods assignment 1
// Max Orok, September 2019
//
//-----------------------------------------------------------------------------
// Godunov's method for linear PDEs (1D)
//-----------------------------------------------------------------------------
// discretize domain into non-overlapping cells
// compute average of initial condition in each cell
// solve a Riemann problem on each boundary to find the fluxes
// advance each cell average
//
//-----------------------------------------------------------------------------
// Dependencies
//-----------------------------------------------------------------------------
// linear algebra
extern crate nalgebra as na;
// floating-point equality
extern crate approx;

// evenly-spaced floats
extern crate itertools_num;
use itertools_num::linspace;

// index operator overloading
extern crate derive_more;
use derive_more::Index;

// numerical integration
extern crate gauss_quad;
use gauss_quad::Midpoint;


// pub fn make_flux_fn<F> (left: f64, right: f64, flux_fn: Fn(f64, f64) -> f64) {
//     match flux_type {
//         LeftCell    => flux_fn(left),
//         RightCell   => flux_fn(right),
//         CellAverage => 0.5 * (flux_fn(left) + flux_fn(right)),
//         Godunov => flux_fn(riemann_soln(left, right)),
//     }
// }

fn riemann_soln(left: f64, right: f64) -> f64 {

    unimplemented!();
}

// the solution vector
#[derive(Debug)]
pub struct Solution1d {
    // cell values
    cells: Vec<f64>,
    // the discretized domain
    domain: Vec<Boundary1d>,
}

impl Solution1d {
    // constructor
    pub fn new(left: f64, right: f64, num_cells: usize) -> Self {

        // initialize the domain
        let mut domain = Vec::new();

        Solution1d {
            cells: vec![0.0; num_cells], // 0 out the solution vector
            domain,
        }
    }

    // initialize the solution vector with some initial condition expression
    pub fn init<F> (&mut self, ic: F)
        where F: Fn(f64) -> f64 {
            // hardcode averaging scheme for simplicity
            let quad = Midpoint::init(self.cells.len());

    }

    // update the solution using some flux function
    pub fn update<F> (&mut self, flux_fn: F)
        where F: Fn(f64) -> f64 {

    }

}

// a boundary -- a cell interface
#[derive(Debug)]
struct Boundary1d {
    coord: f64,
    left_cell: usize,
    right_cell: usize,
}


//    // create `num_cell` equally sized pieces between left and right
//    pub fn create(left: f64, right: f64, num_cells: usize) -> Self {
//        // one more boundary than number of cells
//        let boundary_coords = linspace::<f64>(left, right, num_cells + 1);
//        let mut boundaries = Vec::new();
//
//        for (idx, coord) in boundary_coords.enumerate() {
//            // (-1 + num_cells is for wraparound indexing)
//            let prev_idx = if idx == 0 {num_cells - 1} else {idx - 1};
//            let next_idx = (idx + 1) % (num_cells + 1);
//            boundaries.push(Boundary1D {
//                coord,
//                left_cell: prev_idx,
//                right_cell: next_idx,
//            });
//        }
//
//        PeriodicDomain1d {
//            boundaries,
//        }
//    }
//}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn periodic() {
        let l = 10;
        let domain = PeriodicDomain1d::even(-10., 10., l);
        let left = 0;
        let right = domain.len() - 1;
        // left periodic
        assert!(domain.prev_idx(left) == right);
        // right periodic
        assert!(domain.next_idx(right) == left);
    }

    #[test]
    pub fn domain_size() {
        let l = 100;
        let domain = PeriodicDomain1d::even(-10., 10., l);
        assert!(domain.len() == l + 1);
    }
}
