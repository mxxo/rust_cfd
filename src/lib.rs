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

// a discrete solution cell
#[derive(Debug)]
struct Cell1d {
    pub value: f64,
    pub left: f64,
    pub right: f64,
}

// a boundary | cell interface
#[derive(Debug)]
struct Boundary1d {
    pub coord: f64,
    pub left_cell: usize,
    pub right_cell: usize,
}

// create `num_cell` equally sized pieces between left and right
fn make_domain(left: f64, right: f64, num_cells: usize) -> Vec<Boundary1d> {
    // one more boundary than number of cells
    let boundary_coords = linspace::<f64>(left, right, num_cells + 1);
    let mut boundaries = Vec::new();

    for (idx, coord) in boundary_coords.enumerate() {
        // num_cells is for wraparound indexing of usize values
        let prev_idx = if idx == 0 {num_cells - 1} else { idx - 1 };
        let next_idx = idx % num_cells;
        boundaries.push(Boundary1d {
            coord,
            left_cell: prev_idx,
            right_cell: next_idx,
        });
    }

    boundaries
}

// the solution vector
#[derive(Debug)]
pub struct Solution1d {
    // cell values
    cells: Vec<Cell1d>,
    // the discretized domain
    domain: Vec<Boundary1d>,
}

impl Solution1d {
    // constructor
    pub fn new(left: f64, right: f64, num_cells: usize) -> Self {
        // initialize the domain
        let domain = make_domain(left, right, num_cells);
        // reserve space for the cell vector
        let mut cells = Vec::with_capacity(num_cells);

        for cell_idx in 0..num_cells {
            cells.push(Cell1d {
                value: 0.0,
                left: domain[cell_idx].coord,
                right: domain[cell_idx + 1].coord,
            });
        }

        Solution1d {
            cells,
            domain,
        }
    }

    // initialize the solution vector using an initial condition expression
    pub fn init<F>(&mut self, ic_fn: F)
    where
        F: Fn(f64) -> f64,
    {
        // hardcode averaging scheme for simplicity
        let quad = Midpoint::init(1000);
        // average initial condition in each cell
        for mut cell in self.cells.iter_mut() {
            cell.value = quad.integrate(cell.left, cell.right, &ic_fn) / (cell.right - cell.left);
        }
    }

    // update the solution using some flux function
    pub fn update<F>(&mut self, flux_fn: F)
    where
        F: Fn(f64, f64) -> f64,
    {

    }

    // print out the solution to csv
    pub fn print(&self) {
        println!("position, value");
        for cell in &self.cells {
            println!("{}, {}", 0.5 * (cell.left + cell.right), cell.value);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    //#[test]
    //pub fn periodic() {
    //    let l = 10;
    //    let domain = PeriodicDomain1d::even(-10., 10., l);
    //    let left = 0;
    //    let right = domain.len() - 1;
    //    // left periodic
    //    assert!(domain.prev_idx(left) == right);
    //    // right periodic
    //    assert!(domain.next_idx(right) == left);
    //}

    #[test]
    pub fn domain_size() {
        let l = 100;
        let domain = make_domain(-10., 10., l);
        assert!(domain.len() == l + 1);
    }
}
