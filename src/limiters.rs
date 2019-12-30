//! Slope limiters for assignment 4
//! -- Max Orok December 2019

use crate::euler::EulerCell1d;
use crate::fluxes::EulerCellDelta;
use std::ops::Sub;

/// A set of limiter coefficients for the 1D Euler equations.
pub struct EulerLimit {
    pub density_limit: f64,
    pub momentum_limit: f64,
    pub energy_limit: f64,
}

// could generalize limiters with a Limiter trait
// pub trait Limiter;

/// The *Van Albada* slope limiter.
pub struct VanAlbada;

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

impl VanAlbada {

    fn state_limiter(middle: EulerCell1d, (left, right): (EulerCell1d, EulerCell1d)) -> EulerLimit {

        let backward_diff = middle - left;
        let forward_diff = right - middle;

        EulerLimit {
            density_limit: Self::limit(backward_diff.density, forward_diff.density),
            momentum_limit: Self::limit(backward_diff.momentum, forward_diff.momentum),
            energy_limit: Self::limit(backward_diff.energy, forward_diff.energy),
        }

        //unimplemented!();
    }

    fn limit(backward_difference: f64, forward_difference: f64) -> f64 {
       (backward_difference * forward_difference) * (backward_difference + forward_difference)
           / (backward_difference * backward_difference + forward_difference * forward_difference + 1e-8)
    }

}

#[cfg(test)]

mod tests {
    extern crate rand;
    use rand::Rng;
    use super::*;

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

}
