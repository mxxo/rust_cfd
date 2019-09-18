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

// the solution vector
#[derive(Debug)]
pub struct Solution1d {
    // cell values
    cells: Vec<Cell1d>,
    // the discretized domain
    domain: PeriodicDomain1d,
}

impl Solution1d {
    // constructor
    pub fn new(left: f64, right: f64, num_cells: usize) -> Self {
        // make a problem domain
        let domain = PeriodicDomain1d::even(left, right, num_cells);

        // initialize the solution vector
        let mut cells = Vec::new();
        cells.reserve(num_cells);

        for i in 0..num_cells {
            cells.push(
                Cell1d {value: 0.0, left: i, right: i+1}
            );
        }

        Solution1d {
            cells,
            domain,
        }
    }

//initialize(&self,
}

// a discrete solution cell
#[derive(Debug)]
struct Cell1d {
    pub value: f64,
    pub left: usize,
    pub right: usize,
}



// a 1D periodic domain
#[derive(Debug, Index)]
struct PeriodicDomain1d {
    boundaries: Vec<f64>
}

impl PeriodicDomain1d {
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
        PeriodicDomain1d {
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
