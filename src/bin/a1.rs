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
extern crate nalgebra; 
extern crate nalgebra_lapack; 
use nalgebra::Matrix2;
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

    // find an appropriate timestep using cfl
    let cfl = 1.0;
    let wave_speed = 2.0; // constant in problem

    let time_step = get_time_step(soln.delta_x, wave_speed, cfl);

    // step through time until t_end

    let mut time = 0.0;
    let t_end = 10.;

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

// the second set of pdes
fn pde2(flux_type: &FluxType) -> (Solution1d, Solution1d) {
    let left = -5.0;
    let right = 5.0;
    let num_cells = 100;

    // ------------------------------------------------------------------
    // equation 1
    // ------------------------------------------------------------------
    let mut soln1 = {
        let mut soln = Solution1d::new(left, right, num_cells);
        let ic = |x: f64| if x.abs() < 2.0 { (-x * x).exp() } else { 0. };
        soln.init(ic);

        soln
    };

    // ------------------------------------------------------------------
    // equation 2
    // ------------------------------------------------------------------
    let mut soln2 = {
        let mut soln = Solution1d::new(left, right, num_cells);
        let ic = |x: f64| f64::sin(x * 2.0 * f64::consts::PI / 5.0);
        soln.init(ic);

        soln
    };
    
    // TODO eigensystem stuff here 

    // constant matrix A from the problem description 
    let a = Matrix2::new(0.0, 1.0,
                         0.5, 0.5); 
    
    let eigensystem = Eigen::new(a, true, true).expect("couldn't solve eigensystem"); 
    
    let e_vals = eigensystem.eigenvalues; 
    let r_matrix = eigensystem.eigenvectors.expect("couldn't find eigenvector set"); 
    let r_inverse = eigensystem.left_eigenvectors.expect("couldn't find eigenvector set"); 

    dbg!(e_vals); 
    dbg!(r_matrix); 
    dbg!(r_inverse); 

    // constant wave speeds in the problem
    let wave_speed_1 = 1.0;
    let wave_speed_2 = -0.5;

    // find an appropriate timestep using cfl
    // take the lowest time_step and step through both at once
    let cfl = 1.0;

    let ts1 = get_time_step(soln1.delta_x, wave_speed_1, cfl); 
    let ts2 = get_time_step(soln2.delta_x, wave_speed_2, cfl);

    let time_step = ts1.min(ts2); 
    
    // 
    // TODO: transform to uncoupled equations
    // 

    let mut time = 0.0;
    let t_end = 10.;

    loop {
        if time + time_step > t_end {
            soln1.update(&flux_type, t_end - time, wave_speed_1);
            soln2.update(&flux_type, t_end - time, wave_speed_2);
        } else {
            soln1.update(&flux_type, time_step, wave_speed_1);
            soln2.update(&flux_type, time_step, wave_speed_2);
        }

        time += time_step;
        if time > t_end {
            break;
        }
    }

    // 
    // TODO: transform back to initial equations
    // 

    // return the solutions
    (soln1, soln2)
}

// helper function for timesteps
fn get_time_step(delta_x: f64, wave_speed: f64, cfl: f64) -> f64 {
    cfl * delta_x / wave_speed.abs() // wave speed absolute value to avoid negative timesteps
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
