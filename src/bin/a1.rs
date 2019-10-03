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
use rust_cfd::{FluxType, Solution1d};

// basic math functions
use std::f64;

// eigen decomposition
extern crate nalgebra as na;
extern crate nalgebra_lapack;

use nalgebra_lapack::Eigen;

// filesystem
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

// prints out two solution files
fn main() {
    // choose the flux function here
    let flux_type = FluxType::LeftCell;
    // let flux_type = FluxType::RightCell;
    // let flux_type = FluxType::CellAverage;
    // let flux_type = FluxType::Riemann;

    // first pde
    let soln1 = pde1(&flux_type);
    let (coords, values) = soln1.export();
    print_data("out/question1.csv", "position, value", &coords, &[&values]);

    // second pde
    let (soln2a, soln2b) = pde2(&flux_type);
    // destructure coords and values
    let (coords, a_values) = soln2a.export();
    let (_, b_values) = soln2b.export();

    print_data(
        "out/question2.csv",
        "position, pde1, pde2",
        &coords,
        &[&a_values, &b_values],
    );
}

// the first pde
fn pde1(flux_type: &FluxType) -> Solution1d {
    // discretize domain
    let left = -5.0;
    let right = 5.0;
    let num_cells = 100;

    // make a new, empty solution
    let mut soln = Solution1d::new(left, right, num_cells);

    // update the solution with initial conditions
    let ic = |x: f64| if x.abs() < 2.0 { (-x * x).exp() } else { 0. };
    soln.init(ic);

    // solve the pde
    let cfl = 0.5;
    let wave_speed = 2.0; // constant in problem
    let time = 10.0;

    // return the pde solution
    step_pde(soln, flux_type, wave_speed, time, cfl)
}

// the second set of pdes
fn pde2(flux_type: &FluxType) -> (Solution1d, Solution1d) {
    let left = -5.0;
    let right = 5.0;
    let num_cells = 100;

    // ------------------------------------------------------------------
    // equation 1 initial conditions
    // ------------------------------------------------------------------
    let u0_ic = {
        let mut soln = Solution1d::new(left, right, num_cells);
        let ic = |x: f64| if x.abs() < 2.0 { (-x * x).exp() } else { 0.0 };
        soln.init(ic);

        let (_, values) = soln.export();
        values
    };

    //dbg!(u0_ic);

    // ------------------------------------------------------------------
    // equation 2 initial conditions
    // ------------------------------------------------------------------
    let u1_ic = {
        let mut soln = Solution1d::new(left, right, num_cells);
        let ic = |x: f64| f64::sin(x * 2.0 * f64::consts::PI / 5.0);
        soln.init(ic);

        let (_, values) = soln.export();
        values
    };

    // ------------------------------------------------------------------
    // decouple equations with eigendecomposition
    // ------------------------------------------------------------------

    // constant matrix A from the problem description
    let a = na::Matrix2::new(0.0, 1.0, 0.5, 0.5);

    let eigensystem = Eigen::new(a, true, true).expect("couldn't solve eigensystem");

    let e_vals = eigensystem.eigenvalues;
    let r_matrix = eigensystem
        .eigenvectors
        .expect("couldn't find eigenvector set");

    let r_inverse = r_matrix.try_inverse().expect("couldn't invert R matrix");

    // println!("{}", r_matrix);
    // println!("{}", r_inverse);

    // calculate transformed initial conditions
    // -- would be nice to make this more general
    let transform_soln = |index| -> (Solution1d, f64) {
        let left_evec = r_inverse.row(index);

        let v_ic: Vec<f64> = u0_ic
            .iter()
            .zip(u1_ic.iter())
            .map(|(&u0, &u1)| {
                let soln_piece = na::Matrix1x2::new(u0, u1);
                left_evec.dot(&soln_piece)
            })
            .collect();

        let mut v = Solution1d::new(left, right, num_cells);
        v.fill(&v_ic);

        (v, e_vals[index])
    };

    let (v0, wavespeed_v0) = transform_soln(0);
    let (v1, wavespeed_v1) = transform_soln(1);

    //(v0, v1)

    // solve both equations in v-space
    let time = 10.0;
    let cfl = 0.99;

    // v0
    let solved_v0 = step_pde(v0, flux_type, wavespeed_v0, time, cfl);
    let (_, v0_values) = solved_v0.export(); 

    // v1
    let solved_v1 = step_pde(v1, flux_type, wavespeed_v1, time, cfl);
    let (_, v1_values) = solved_v1.export(); 

    //(solved_v0, solved_v1)

    // transform back to initial equations
    let transform_back = |index| -> Solution1d { 
        let e_vec = r_matrix.row(index); 
        let v_ic: Vec<f64> = v0_values
            .iter()
            .zip(v1_values.iter())
            .map(|(&v0, &v1)| {
                let soln_piece = na::Matrix1x2::new(v0, v1);
                e_vec.dot(&soln_piece)
            })
            .collect();

        let mut v = Solution1d::new(left, right, num_cells);
        v.fill(&v_ic);

        v
    };

    let u0 = transform_back(0);
    let u1 = transform_back(1);

    // return the solutions
    (u0, u1)
}

// helper method to solve pdes
// -- method is kinda smelly... 3 floats in a row
fn step_pde(
    mut soln: Solution1d,
    flux_type: &FluxType,
    wave_speed: f64,
    t_end: f64,
    cfl: f64,
) -> Solution1d {
    // find appropriate timestep for this solution
    let time_step = get_time_step(soln.delta_x, wave_speed, cfl);

    // step through time until t_end
    let mut time = 0.0;

    loop {
        if time + time_step > t_end {
            soln.update(&flux_type, t_end - time, wave_speed);
        } else {
            soln.update(&flux_type, time_step, wave_speed);
        }

        time += time_step;
        if time > t_end {
            break;
        }
    }

    soln
}

// helper function for timesteps
fn get_time_step(delta_x: f64, wave_speed: f64, cfl: f64) -> f64 {
    // wave speed absolute value to avoid negative timesteps
    cfl * delta_x / wave_speed.abs()
}

// simple csv printing function
fn print_data(filename: &str, header: &str, x_values: &[f64], y_values: &[&[f64]]) {
    // unwraps are quick and dirty for now

    let path = Path::new(filename);
    let mut filestream = BufWriter::new(File::create(&path).unwrap());

    write!(&mut filestream, "{}", header).unwrap();

    for (index, x_val) in x_values.iter().enumerate() {
        write!(&mut filestream, "\n{}", x_val).unwrap();
        for y_column in y_values {
            write!(&mut filestream, ", {}", y_column[index]).unwrap();
        }
    }
}
