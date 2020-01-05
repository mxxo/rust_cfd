//! 2d approximate fluxes for the Euler equations.

use super::*;
use crate::euler2d::*;

pub trait FluxFunction2d: Copy {
    fn calculate_x_flux(&self, left: EulerCell2d, right: EulerCell2d) -> EulerFlux2d;
    fn calculate_y_flux(&self, bottom: EulerCell2d, top: EulerCell2d) -> EulerFlux2d;
}

#[derive(Debug, Clone, Copy)]
pub struct Hlle2d;

impl Hlle2d {

    pub fn middle_x_flux(
        left: StateVec2d,
        right: StateVec2d,
        l_minus: f64,
        l_plus: f64,
    ) -> EulerFlux2d {
        let lr_term =
            1.0 / (l_plus - l_minus) * (l_plus * left.to_x_flux() - l_minus * right.to_x_flux());

        let u_delta = right - left;
        let delta_term = (l_plus * l_minus) / (l_plus - l_minus)
            * EulerFlux2d {
                /* dummy EulerFlux from cell deltas */
                density_flux: u_delta.density,
                x_momentum_flux: u_delta.x_momentum,
                y_momentum_flux: u_delta.y_momentum,
                energy_flux: u_delta.energy,
            };

        lr_term + delta_term
    }

    pub fn middle_y_flux(
        bottom: StateVec2d,
        top: StateVec2d,
        l_minus: f64,
        l_plus: f64,
    ) -> EulerFlux2d {
        let lr_term =
            1.0 / (l_plus - l_minus) * (l_plus * bottom.to_y_flux() - l_minus * top.to_y_flux());

        let u_delta = top - bottom;
        let delta_term = (l_plus * l_minus) / (l_plus - l_minus)
            * EulerFlux2d {
                /* dummy EulerFlux from cell deltas */
                density_flux: u_delta.density,
                x_momentum_flux: u_delta.x_momentum,
                y_momentum_flux: u_delta.y_momentum,
                energy_flux: u_delta.energy,
            };

        lr_term + delta_term
    }


    // x-fluxes

    pub fn left_wavespeed(left: EulerPrimitive2d, roe_avg: EulerPrimitive2d) -> f64 {
        f64::min(
            left.x_vel - left.sound_speed(),
            roe_avg.x_vel - roe_avg.sound_speed(),
        )
    }

    pub fn right_wavespeed(right: EulerPrimitive2d, roe_avg: EulerPrimitive2d) -> f64 {
        f64::max(
            right.x_vel + right.sound_speed(),
            roe_avg.x_vel + roe_avg.sound_speed(),
        )
    }

    // y-fluxes

    pub fn bottom_wavespeed(bottom: EulerPrimitive2d, roe_avg: EulerPrimitive2d) -> f64 {
        f64::min(
            bottom.y_vel - bottom.sound_speed(),
            roe_avg.y_vel - roe_avg.sound_speed(),
        )
    }

    pub fn top_wavespeed(top: EulerPrimitive2d, roe_avg: EulerPrimitive2d) -> f64 {
        f64::max(
            top.y_vel + top.sound_speed(),
            roe_avg.y_vel + roe_avg.sound_speed(),
        )
    }
}

impl FluxFunction2d for Hlle2d {
    fn calculate_x_flux(&self, left: EulerCell2d, right: EulerCell2d) -> EulerFlux2d {

        // find roe average along x-direction
        let roe_avg = Roe2d::average_state(left.state, right.state);
        let l_minus = Hlle2d::left_wavespeed(left.state.into(), roe_avg);
        let l_plus = Hlle2d::right_wavespeed(right.state.into(), roe_avg);

        if l_minus > 0.0 {
            left.state.to_x_flux()
        } else if l_plus < 0.0 {
            right.state.to_x_flux()
        } else {
            Hlle2d::middle_x_flux(left.state, right.state, l_minus, l_plus)
        }
    }

    fn calculate_y_flux(&self, bottom: EulerCell2d, top: EulerCell2d) -> EulerFlux2d {

        // find roe average along y-direction
        let roe_avg = Roe2d::average_state(bottom.state, top.state);
        let l_minus = Hlle2d::bottom_wavespeed(bottom.state.into(), roe_avg.into());
        let l_plus = Hlle2d::top_wavespeed(top.state.into(), roe_avg.into());

        if l_minus > 0.0 {
            bottom.state.to_y_flux()
        } else if l_plus < 0.0 {
            top.state.to_y_flux()
        } else {
            Hlle2d::middle_y_flux(bottom.state, top.state, l_minus, l_plus)
        }
    }
}



/// Two dimensional Roe equations.
#[derive(Debug, Clone, Copy)]
pub struct Roe2d;

// quick and dirty reuse of the 1d functions
impl EulerPrimitive2d {
    /// Get a view of this 2d state as a 1d state in the x-direction.
    pub fn as_x_primitive(self) -> PrimitiveState {
        PrimitiveState {
            density: self.density,
            velocity: self.x_vel,
            pressure: self.pressure,
            gamma: self.gamma,
        }
    }

    /// Get a view of this 2d state as a 1d state in the y-direction.
    pub fn as_y_primitive(self) -> PrimitiveState {
        PrimitiveState {
            density: self.density,
            velocity: self.y_vel,
            pressure: self.pressure,
            gamma: self.gamma,
        }
    }

    /// Get a view of this state using the combined velocities.
    pub fn as_xy_primitive(self) -> PrimitiveState {
        PrimitiveState {
            density: self.density,
            velocity: self.velocity_sq().sqrt(),
            pressure: self.pressure,
            gamma: self.gamma,
        }
    }
}

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

#[cfg(test)]
mod tests {
    use super::*;
}
