//! Slope limiters for assignment 4
//! -- Max Orok December 2019

use crate::euler::EulerCell1d;
use crate::fluxes::EulerCellDelta;
use std::ops::{Mul, Sub};

/// A set of limiter coefficients for the 1D Euler equations.

#[derive(Debug, Copy, Clone)]
pub struct EulerLimit {
    pub density_limit: f64,
    pub momentum_limit: f64,
    pub energy_limit: f64,
}

// could generalize limiters with a Limiter trait
// pub trait Limiter;

/// The *Van Albada* slope limiter.
#[derive(Debug, Copy, Clone)]
pub struct VanAlbada;

impl EulerCell1d {
    /// Reconstruct the state using a limiter. Scale the limiter by -1.0 as required.
    /// E.g. u_left at the left interface is u_i-1 + limiter_i-1
    /// but u_right at the left interface is u_i - limiter_i.
    pub fn reconstruct(self, limiter: EulerLimit) -> EulerCell1d {
        let limiter = 0.5 * self.width() * limiter;

        EulerCell1d {
            density: self.density + limiter.density_limit,
            momentum: self.momentum + limiter.momentum_limit,
            energy: self.energy + limiter.energy_limit,
            ..self
        }
    }
}

/// A ad-hoc way to get EulerCellDelta values from EulerCells.
impl Sub for EulerCell1d {
    type Output = EulerCellDelta;

    fn sub(self, other: Self) -> Self::Output {
        Self::Output {
            density: self.density - other.density,
            momentum: self.momentum - other.momentum,
            energy: self.energy - other.energy,
        }
    }
}

impl Mul<EulerCellDelta> for f64 {
    type Output = EulerCellDelta;

    fn mul(self, rhs: Self::Output) -> Self::Output {
        Self::Output {
            density: self * rhs.density,
            momentum: self * rhs.momentum,
            energy: self * rhs.energy,
        }
    }
}

impl Mul<EulerLimit> for f64 {
    type Output = EulerLimit;

    fn mul(self, rhs: Self::Output) -> Self::Output {
        Self::Output {
            density_limit: self * rhs.density_limit,
            momentum_limit: self * rhs.momentum_limit,
            energy_limit: self * rhs.energy_limit,
        }
    }
}

impl VanAlbada {
    /// Get the series of slope limiters for a solution.
    pub fn soln_limiters(soln: &Vec<EulerCell1d>) -> Vec<EulerLimit> {
        let mut limiters: Vec<EulerLimit> = Vec::with_capacity(soln.len());

        // first cell take a copy as left cell
        limiters.push(Self::state_limiter(soln[0], (soln[0], soln[1])));

        for i in 1..soln.len() - 1 {
            limiters.push(Self::state_limiter(soln[i], (soln[i - 1], soln[i + 1])));
        }

        // last cell take a copy as right cell
        let last_idx = soln.len() - 1;
        limiters.push(Self::state_limiter(
            soln[last_idx],
            (soln[last_idx - 1], soln[last_idx]),
        ));

        limiters
    }

    fn state_limiter(middle: EulerCell1d, (left, right): (EulerCell1d, EulerCell1d)) -> EulerLimit {
        let backward_diff = 1.0 / middle.width() * (middle - left);
        let forward_diff = 1.0 / middle.width() * (right - middle);

        EulerLimit {
            density_limit: Self::limit(backward_diff.density, forward_diff.density),
            momentum_limit: Self::limit(backward_diff.momentum, forward_diff.momentum),
            energy_limit: Self::limit(backward_diff.energy, forward_diff.energy),
        }
    }

    #[inline(always)]
    fn limit(backward_difference: f64, forward_difference: f64) -> f64 {
        (backward_difference * forward_difference) * (backward_difference + forward_difference)
            / (backward_difference * backward_difference
                + forward_difference * forward_difference
                + 1e-8)
    }
}

#[cfg(test)]

mod tests {
    extern crate approx;
    extern crate rand;

    use super::*;
    use crate::euler::PrimitiveState;
    use approx::*;
    use rand::Rng;

    #[test]
    fn reconstruct_bounds() {
        let soln = vec![
            EulerCell1d::new(
                /*let left_state = */
                PrimitiveState {
                    density: 2.0,
                    velocity: 164.83,
                    pressure: 201.17e3,
                    gamma: 1.4,
                },
                (0.0, 1.0),
            ),
            EulerCell1d::new(
                /*let left_state = */
                PrimitiveState {
                    density: 40.0,
                    velocity: 200.0,
                    pressure: 400e3,
                    gamma: 1.4,
                },
                (1.0, 2.0),
            ),
            EulerCell1d::new(
                /*let left_state = */
                PrimitiveState {
                    density: 60.0,
                    velocity: 300.0,
                    pressure: 500e3,
                    gamma: 1.4,
                },
                (2.0, 3.0),
            ),
        ];

        let limiters = VanAlbada::soln_limiters(&soln);

        let left_cell = soln[0];
        let middle_cell = soln[1];
        let right_cell = soln[2];

        // -- left interface
        let left_reconstruct = left_cell.reconstruct(limiters[0]);
        let right_reconstruct = middle_cell.reconstruct(-1.0 * limiters[1]);

        // - left and left reconstruct should be the same
        assert_relative_eq!(left_cell.density, left_reconstruct.density);
        assert_relative_eq!(left_cell.momentum, left_reconstruct.momentum);
        assert_relative_eq!(left_cell.energy, left_reconstruct.energy);

        // - right reconstruct should be larger than left_cell but smaller than middle cell
        assert!(right_reconstruct.density > left_cell.density);
        assert!(right_reconstruct.density < middle_cell.density);
        assert!(right_reconstruct.momentum > left_cell.momentum);
        assert!(right_reconstruct.momentum < middle_cell.momentum);
        assert!(right_reconstruct.energy > left_cell.energy);
        assert!(right_reconstruct.energy < middle_cell.energy);

        // -- right interface
        let left_reconstruct = middle_cell.reconstruct(limiters[1]);
        let right_reconstruct = right_cell.reconstruct(-1.0 * limiters[2]);

        // - left reconstruct should be larger than middle_cell but smaller than right cell
        assert!(left_reconstruct.density > middle_cell.density);
        assert!(left_reconstruct.density < right_cell.density);
        assert!(left_reconstruct.momentum > middle_cell.momentum);
        assert!(left_reconstruct.momentum < right_cell.momentum);
        assert!(left_reconstruct.energy > middle_cell.energy);
        assert!(left_reconstruct.energy < right_cell.energy);

        // - right and right reconstruct should be the same
        assert_relative_eq!(right_cell.density, right_reconstruct.density);
        assert_relative_eq!(right_cell.momentum, right_reconstruct.momentum);
        assert_relative_eq!(right_cell.energy, right_reconstruct.energy);
    }

    #[test]
    fn reconstruct_zero() {
        let soln = vec![
            EulerCell1d::new(
                /*let left_state = */
                PrimitiveState {
                    density: 2.281,
                    velocity: 164.83,
                    pressure: 201.17e3,
                    gamma: 1.4,
                },
                (0.0, 1.0),
            ),
            EulerCell1d::new(
                /*let left_state = */
                PrimitiveState {
                    density: 2.281,
                    velocity: 164.83,
                    pressure: 201.17e3,
                    gamma: 1.4,
                },
                (1.0, 2.0),
            ),
            EulerCell1d::new(
                /*let left_state = */
                PrimitiveState {
                    density: 2.281,
                    velocity: 164.83,
                    pressure: 201.17e3,
                    gamma: 1.4,
                },
                (2.0, 3.0),
            ),
        ];

        let limiters = VanAlbada::soln_limiters(&soln);

        let left_cell = soln[0];
        let middle_cell = soln[1];
        let right_cell = soln[2];

        let left_reconstruct = left_cell.reconstruct(limiters[0]);
        let right_reconstruct = middle_cell.reconstruct(-1.0 * limiters[1]);

        dbg!(&left_reconstruct);
        dbg!(&right_reconstruct);

        // left interface
        // - left and left reconstruct should be the same
        // - middle and right reconstruct should be the same
        assert_relative_eq!(left_cell.density, left_reconstruct.density);
        assert_relative_eq!(left_cell.momentum, left_reconstruct.momentum);
        assert_relative_eq!(left_cell.energy, left_reconstruct.energy);

        assert_relative_eq!(middle_cell.density, right_reconstruct.density);
        assert_relative_eq!(middle_cell.momentum, right_reconstruct.momentum);
        assert_relative_eq!(middle_cell.energy, right_reconstruct.energy);
    }

    #[test]
    fn probably_tvd() {
        let mut rng = rand::thread_rng();
        for _ in 0..100 {
            let limit = VanAlbada::limit(rng.gen::<f64>(), rng.gen::<f64>());
            assert!(limit <= 1.0 && limit >= 0.0);
        }
    }

    #[test]
    fn flat_limiters() {
        let mut rng = rand::thread_rng();

        let limit = VanAlbada::limit(0.0, 1.0);
        assert!(limit < 1e-6);
        let limit = VanAlbada::limit(1.0, 0.0);
        assert!(limit < 1e-6);

        for _ in 0..100 {
            let rando = rng.gen::<f64>();
            let limit = VanAlbada::limit(rando, 0.0);
            assert!(limit < 1e-6);
            let limit = VanAlbada::limit(0.0, rando);
            assert!(limit < 1e-6);
        }
    }

    #[test]
    fn zero_limiters() {
        let soln = vec![
            EulerCell1d::new(
                /*let left_state = */
                PrimitiveState {
                    density: 2.281,
                    velocity: 164.83,
                    pressure: 201.17e3,
                    gamma: 1.4,
                },
                (0.0, 1.0),
            ),
            EulerCell1d::new(
                /*let left_state = */
                PrimitiveState {
                    density: 2.281,
                    velocity: 164.83,
                    pressure: 201.17e3,
                    gamma: 1.4,
                },
                (1.0, 2.0),
            ),
            EulerCell1d::new(
                /*let left_state = */
                PrimitiveState {
                    density: 2.281,
                    velocity: 164.83,
                    pressure: 201.17e3,
                    gamma: 1.4,
                },
                (2.0, 3.0),
            ),
        ];

        let limiters = VanAlbada::soln_limiters(&soln);

        assert_eq!(soln.len(), limiters.len());

        for limiter in limiters {
            assert_relative_eq!(limiter.density_limit, 0.0);
            assert_relative_eq!(limiter.momentum_limit, 0.0);
            assert_relative_eq!(limiter.energy_limit, 0.0);
        }
    }

    #[test]
    // Higher state in the middle, all limiters should be zero to avoid introducing
    // new extrema.
    fn peak_limiters() {
        let soln = vec![
            EulerCell1d::new(
                /*let left_state = */
                PrimitiveState {
                    density: 2.281,
                    velocity: 164.83,
                    pressure: 201.17e3,
                    gamma: 1.4,
                },
                (0.0, 1.0),
            ),
            EulerCell1d::new(
                /*let left_state = */
                PrimitiveState {
                    density: 4.0,
                    velocity: 200.0,
                    pressure: 400e3,
                    gamma: 1.4,
                },
                (1.0, 2.0),
            ),
            EulerCell1d::new(
                /*let left_state = */
                PrimitiveState {
                    density: 2.281,
                    velocity: 164.83,
                    pressure: 201.17e3,
                    gamma: 1.4,
                },
                (2.0, 3.0),
            ),
        ];

        let limiters = VanAlbada::soln_limiters(&soln);

        assert_eq!(soln.len(), limiters.len());

        for limiter in limiters {
            assert_relative_eq!(limiter.density_limit, 0.0);
            assert_relative_eq!(limiter.momentum_limit, 0.0);
            assert_relative_eq!(limiter.energy_limit, 0.0);
        }
    }

    #[test]
    // States increase from left to right, limiters should be nonzero.
    fn step_limiters() {
        let soln = vec![
            EulerCell1d::new(
                /*let left_state = */
                PrimitiveState {
                    density: 2.0,
                    velocity: 164.83,
                    pressure: 201.17e3,
                    gamma: 1.4,
                },
                (0.0, 1.0),
            ),
            EulerCell1d::new(
                /*let left_state = */
                PrimitiveState {
                    density: 40.0,
                    velocity: 200.0,
                    pressure: 400e3,
                    gamma: 1.4,
                },
                (1.0, 2.0),
            ),
            EulerCell1d::new(
                /*let left_state = */
                PrimitiveState {
                    density: 60.0,
                    velocity: 300.0,
                    pressure: 500e3,
                    gamma: 1.4,
                },
                (2.0, 3.0),
            ),
        ];

        let limiters = VanAlbada::soln_limiters(&soln);

        assert_eq!(soln.len(), limiters.len());

        dbg!(&limiters);
        dbg!(&limiters[1..2]);

        for limiter in &limiters[1..2] {
            assert!(limiter.density_limit > 10.0);
            assert!(limiter.momentum_limit > 10.0);
            assert!(limiter.energy_limit > 10.0);
        }
    }
}
