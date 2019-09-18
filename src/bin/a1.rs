// Finite-Volume methods assignment 1
// Max Orok, September 2019
//-----------------------------------------------------------------------------
// Godunov's method for linear PDEs (1D)
//-----------------------------------------------------------------------------
// discretize domain into non-overlapping cells
// compute average of initial condition in each cell
// solve a Riemann problem on each boundary to find the fluxes
// advance each cell average

// general solution code
use rust_cfd::{Solution1d};

// e
use std::f64;

extern crate gauss_quad;
use gauss_quad::Midpoint;

// specific solution
fn main() {

    // discretize domain
    let left = -5.0;
    let right = 5.0;
    let num_cells = 3;

    // make a new, empty solution
    let soln = Solution1d::new(left, right, num_cells);
    dbg!(soln);

    // PDE number 1
    // update the solution with initial conditions
    let ic_1 = |x: f64| if x.abs() < 2.0 {
                            (-x * x).exp()
                        } else {
                            0.
                        };

    dbg!(ic_1(2.));
    dbg!(ic_1(0.2));
    dbg!(ic_1(0.));

    let quad = Midpoint::init(10);
    dbg!(quad.integrate(-1.0, 1.0, |x| x*x));

    //soln.update();

    let t_end = 10.;
}
