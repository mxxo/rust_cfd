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

// the solution vector
#[derive(Debug)]
pub struct Solution1d {
    // cell values
    cells: Vec<Cell1d>,
    // the discretized domain
    domain: Domain1d,
}

impl Solution1d {
    // constructor
    pub fn new(left: f64, right: f64, num_cells: usize) -> Self {

        // make a problem domain
        let d = Domain1d::even(left, right, num_cells);
        Solution1d {
            cells: Vec::new(),
            domain: d,
        }
    }
}

// a discrete solution cell
#[derive(Debug)]
struct Cell1d {
    pub value: f64,
    pub left: usize,
    pub right: usize,
}

// a 1D domain
#[derive(Debug, Index)]
struct Domain1d {
    boundaries: Vec<f64>
}

impl Domain1d {
    // number of discrete points
    pub fn len(&self) -> usize {
        self.boundaries.len()
    }

    pub fn next_idx(&self, index: usize) -> usize {
        (index + 1) % self.len()
    }

    pub fn prev_idx(&self, index: usize) -> usize {
        if index == 0 {
            self.len() - 1
        } else {
            index - 1
        }
    }

    // create `num_cell` equally sized pieces between left and right
    pub fn even(left: f64, right: f64, num_cells: usize) -> Self {
        // one more boundary than number of cells
        let boundaries = linspace::<f64>(left, right, num_cells + 1).collect();
        Domain1d {
            boundaries,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn periodic() {
        let l = 10;
        let domain = Domain1d::even(-10., 10., l);
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
        let domain = Domain1d::even(-10., 10., l);
        assert!(domain.len() == l + 1);
    }
}
