//! 2d approximate fluxes for the Euler equations.

use crate::euler2d::*;
use super::*;

#[derive(Debug, Clone, Copy)]
pub struct Hlle2d;

pub trait FluxFunction2d {
    fn calculate_x_flux(&self, left: EulerCell2d, right: EulerCell2d, time_step: f64) -> EulerFlux2d;
    fn calculate_y_flux(&self, bottom: EulerCell2d, top: EulerCell2d, time_step: f64) -> EulerFlux2d;
}

impl FluxFunction2d for Hlle2d {
    fn calculate_x_flux(&self, left: EulerCell2d, right: EulerCell2d, time_step: f64) -> EulerFlux2d {
        unimplemented!();
    }

    fn calculate_y_flux(&self, bottom: EulerCell2d, top: EulerCell2d, time_step: f64) -> EulerFlux2d {
        unimplemented!();
    }
}

/// Two dimensional Roe equations.
#[derive(Debug, Clone, Copy)]
pub struct Roe2d;

// 2d helper functions for Roe average state
impl Roe2d {

    pub fn average_state<P>(left: P, right: P) -> EulerPrimitive2d
    where
        P: Into<EulerPrimitive2d>,
    {
        let left: EulerPrimitive2d = left.into();
        let right: EulerPrimitive2d = right.into();

        EulerPrimitive2d {
            density: Roe::average_density(left.as_x_primitive(), right.as_x_primitive()),
            x_vel: Roe::average_velocity(left.as_x_primitive(), right.as_x_primitive()),
            y_vel: Roe::average_velocity(left.as_y_primitive(), right.as_y_primitive()),
            pressure: Roe::average_pressure(left.as_xy_primitive(), right.as_xy_primitive()),
            gamma: left.gamma,
        }
    }
}


