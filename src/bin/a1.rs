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

// specific solution
fn main() {

    // discretize domain
    let left = -5.0;
    let right = 5.0;
    let num_cells = 10;

    // make a new, empty solution
    let soln_0 = Solution1d::new(left, right, num_cells);
    dbg!(soln_0);

    // set the initial conditions
    //let initial_conditions =

    let t_end = 10.;
}
