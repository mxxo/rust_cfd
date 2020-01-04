//! Assignment 5 Max Orok
//! December-January 2019
//!
//! -----------------------------------------------------------------------------
//! Approximate 2D Euler solver
//! -----------------------------------------------------------------------------

// use rust_cfd::euler1d::{EulerSolution1d, PrimitiveResult, PrimitiveState};
use fluxes::FluxFunction;
use rust_cfd::euler2d::*;
use rust_cfd::fluxes;

// constant for all gases
const GAMMA: f64 = 1.4;
const MIN: Point2d = Point2d { x: -0.5, y: -0.5 };
const MAX: Point2d = Point2d { x: 0.5, y: 0.5 };

fn main() {
    let sq_width = 500;
    let cfl = 0.5;
    let flux_fn = fluxes::Hlle {};

    case_1(sq_width, flux_fn, cfl);
    case_2(sq_width, flux_fn, cfl);
    case_3(sq_width, flux_fn, cfl);
}

fn case_1(sq_width: usize, flux_fn: impl FluxFunction, cfl: f64) {
    let mut soln = EulerSolution2d::square(sq_width, MIN, MAX);

    let u1 = EulerPrimitive2d {
        density: 1.225,
        x_vel: 0.0,
        y_vel: 0.0,
        pressure: 101_325.0,
        gamma: GAMMA,
    };

    let u2 = EulerPrimitive2d {
        density: 0.30625,
        x_vel: 0.0,
        y_vel: 0.0,
        pressure: 25_331.25,
        gamma: GAMMA,
    };

    let case_1_init = |point: Point2d| {
        if point.x < 0.0 && point.y < 0.0 {
            u2
        } else {
            u1
        }
    };

    // initialize solution grid
    soln.init(case_1_init);

    let t_final = 0.75e-3; // 0.75 ms

    soln.write_gmsh("data/case1.msh");
}

fn case_2(sq_width: usize, flux_fn: impl FluxFunction, cfl: f64) {
    let mut soln = EulerSolution2d::square(sq_width, MIN, MAX);

    let u1 = EulerPrimitive2d {
        density: 1.5,
        x_vel: 0.0,
        y_vel: 0.0,
        pressure: 150_000.0,
        gamma: GAMMA,
    };

    let u2 = EulerPrimitive2d {
        density: 0.5323,
        x_vel: 381.3850,
        y_vel: 0.0,
        pressure: 30_000.0,
        gamma: GAMMA,
    };

    let u3 = EulerPrimitive2d {
        density: 0.1380,
        x_vel: 381.3850,
        y_vel: 381.3850,
        pressure: 2_903.2,
        gamma: GAMMA,
    };

    let u4 = EulerPrimitive2d {
        density: 0.5323,
        x_vel: 0.0,
        y_vel: 381.3850,
        pressure: 30_000.0,
        gamma: GAMMA,
    };

    let init_fn = |point: Point2d| {
        if point.x > 0.3 && point.y > 0.3 {
            u1
        } else if point.x < 0.3 && point.y > 0.3 {
            u2
        } else if point.x < 0.3 && point.y < 0.3 {
            u3
        } else {
            u4
        }
    };

    // initialize solution grid
    soln.init(init_fn);

    let t_final = 2.53e-3; // s

    soln.write_gmsh("data/case2.msh");
}

fn case_3(sq_width: usize, flux_fn: impl FluxFunction, cfl: f64) {
    let mut soln = EulerSolution2d::square(sq_width, MIN, MAX);

    let u1 = EulerPrimitive2d {
        density: 1.0,
        x_vel: 237.171,
        y_vel: -158.114,
        pressure: 100_000.0,
        gamma: GAMMA,
    };

    let u2 = EulerPrimitive2d {
        density: 3.0,
        x_vel: 237.171,
        y_vel: 158.114,
        pressure: 100_000.0,
        gamma: GAMMA,
    };

    let u3 = EulerPrimitive2d {
        density: 1.0,
        x_vel: -237.171,
        y_vel: 158.114,
        pressure: 100_000.0,
        gamma: GAMMA,
    };

    let u4 = EulerPrimitive2d {
        density: 3.0,
        x_vel: -237.171,
        y_vel: -158.114,
        pressure: 100_000.0,
        gamma: GAMMA,
    };

    let init_fn = |point: Point2d| {
        if point.x > 0.0 && point.y > 0.0 {
            u1
        } else if point.x < 0.0 && point.y > 0.0 {
            u2
        } else if point.x < 0.0 && point.y < 0.0 {
            u3
        } else {
            u4
        }
    };

    // initialize solution grid
    soln.init(init_fn);

    let t_final = 1.9e-3; // s

    soln.write_gmsh("data/case3.msh");
}
