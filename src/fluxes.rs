/// Euler flux functions for assignment 3

/// The concept of a flux is something that takes two cells and finds the flux
/// at the boundary.
use std::ops::Deref;

use crate::euler::{EulerCell1d, EulerFlux, PrimitiveState};
use crate::riemann::{solve_euler, DomainBounds, StateSide};

pub trait FluxFunction {
    fn calculate_flux(&self, left: EulerCell1d, right: EulerCell1d, time_step: f64) -> EulerFlux;
}

// assignment 3
/// The four fluxes we'll examine for assignment 3.
pub mod first_order {

    use super::*;

    #[derive(Debug, Clone, Copy)]
    /// The flux using an exact solution to the Euler equations.
    pub struct Exact;

    impl FluxFunction for Exact {
        // solve the exact riemann problem
        fn calculate_flux(
            &self,
            left_cell: EulerCell1d,
            right_cell: EulerCell1d,
            time_step: f64,
        ) -> EulerFlux {
            // pose initial conditions based on cells
            let left_state = left_cell.to_euler_state(StateSide::Left);
            let right_state = right_cell.to_euler_state(StateSide::Right);

            let soln = solve_euler(left_state, right_state);

            // dbg!(&soln);

            let domain = DomainBounds {
                left: left_cell.left,
                interface: left_cell.right,
                right: right_cell.right,
            };

            let point_soln = soln.reconstruct_at_point(
                left_state,
                right_state,
                domain,
                time_step,
                domain.interface,
            );

            // todo fix for different gammas
            EulerFlux::new(PrimitiveState::from_DataPoint(point_soln, left_state.gamma))
        }
    }

    #[derive(Debug, Clone, Copy)]
    pub struct Hlle;

    #[derive(Debug, Clone, Copy)]
    pub struct Roe;

    impl PrimitiveState {
        pub fn enthalpy(&self) -> f64 {
            self.gamma * self.pressure / (self.density * (self.gamma - 1.0))
                + (self.velocity * self.velocity / 2.0)
        }
    }

    impl Roe {

        pub fn average_state(left: PrimitiveState, right: PrimitiveState) -> PrimitiveState {
            PrimitiveState {
                density: Self::average_density(left, right),
                velocity: Self::average_velocity(left, right),
                pressure: Self::average_pressure(left, right),
                gamma: left.gamma,
            }
        }

        pub fn average_density(left: PrimitiveState, right: PrimitiveState) -> f64 {
            (left.density * right.density).sqrt()
        }

        pub fn average_velocity(left: PrimitiveState, right: PrimitiveState) -> f64 {
            (left.density.sqrt() * left.velocity + right.density.sqrt() * right.velocity)
                / (left.density.sqrt() + right.density.sqrt())
        }

        pub fn average_pressure(left: PrimitiveState, right: PrimitiveState) -> f64 {
            let avg_enth = Self::average_enthalpy(left, right);
            let avg_vel = Self::average_velocity(left, right);

            (avg_enth - (avg_vel * avg_vel / 2.0))
                * (Self::average_density(left, right) * (left.gamma - 1.0))
                / left.gamma
        }

        pub fn average_enthalpy(left: PrimitiveState, right: PrimitiveState) -> f64 {
            (left.density.sqrt() * left.enthalpy() + right.density.sqrt() * right.enthalpy())
                / (left.density.sqrt() + right.density.sqrt())
        }


        /// Inner flux is -0.5 * sum(a_i, | /\_i |, K_i)
        fn inner_flux(left: EulerCell1d, right: EulerCell1d) -> EulerFlux {

             let state_delta = EulerCellDelta {
                density: right.density - left.density,
                momentum: right.momentum - left.momentum,
                energy: right.energy - left.energy,
            };

            let roe_avg_state = Roe::average_state(left.to_primitive(), right.to_primitive());

            let eigenvalues = (roe_avg_state.velocity - roe_avg_state.sound_speed(),
                               roe_avg_state.velocity,
                               roe_avg_state.velocity + roe_avg_state.sound_speed(),
                               );

            dbg!(eigenvalues);

            unimplemented!();
        }

    }


   struct EulerCellDelta {
       pub density: f64,
       pub momentum: f64,
       pub energy: f64,
   }

    impl FluxFunction for Roe {

        /// Note: this implementation follows Toro's algorithm (Riemann Solvers, p. 353).
        fn calculate_flux(&self, left: EulerCell1d, right: EulerCell1d, _time_step: f64) -> EulerFlux {
            left.to_flux() + right.to_flux() + Roe::inner_flux(left, right)
        }

    }
    // Roe with entropy fix is implemented as a thin wrapper around Roe
    #[derive(Debug, Clone, Copy)]
    pub struct RoeEntropyFix(Roe);


    impl Deref for RoeEntropyFix {
        type Target = Roe;

        fn deref(&self) -> &Roe {
            &self.0
        }
    }

    #[cfg(test)]
    mod tests {
        use approx::*;
        use super::*;

        #[test]
        fn roe_average() {

            let left_state = PrimitiveState {
                density: 2.281,
                velocity: 164.83,
                pressure: 201.17e3,
                gamma: 1.4,
            };

            assert_relative_eq!(Roe::average_density(left_state, left_state.clone()), left_state.density, epsilon = 1e-6);
            assert_relative_eq!(Roe::average_velocity(left_state, left_state.clone()), left_state.velocity, epsilon = 1e-6);
            assert_relative_eq!(Roe::average_pressure(left_state, left_state.clone()), left_state.pressure, epsilon = 1e-6);
            assert_relative_eq!(Roe::average_enthalpy(left_state, left_state.clone()), left_state.enthalpy(), epsilon = 1e-6);
        }

    }

}

pub mod second_order {
    use super::*;
}
