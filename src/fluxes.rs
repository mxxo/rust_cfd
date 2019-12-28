//! Euler flux functions for assignment 3
//!
//! The concept of a flux is something that takes two cells and finds the flux
//! at the boundary.

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

    impl FluxFunction for Roe {
        /// Note: this implementation follows Toro's algorithm (Riemann Solvers, p. 353).
        fn calculate_flux(
            &self,
            left: EulerCell1d,
            right: EulerCell1d,
            _time_step: f64,
        ) -> EulerFlux {
           let roe_avg = Roe::average_state(left.to_primitive(), right.to_primitive());

            0.5 * (left.to_flux()
                    + right.to_flux())
                + Roe::inner_flux(
                    roe_avg,
                    right.state_delta(left),
                    Roe::eigenvalues(roe_avg),
                )
        }
    }

    impl PrimitiveState {
        pub fn enthalpy(&self) -> f64 {
            self.gamma * self.pressure / (self.density * (self.gamma - 1.0))
                + (self.velocity * self.velocity / 2.0)
        }

        /// Returns (left, center, right) state eigenvalues.
        pub fn eigenvalues(&self) -> (f64, f64, f64) {
            (
                self.velocity - self.sound_speed(),
                self.velocity,
                self.velocity + self.sound_speed(),
            )
        }
    }

    impl EulerCell1d {
        pub(crate) fn state_delta(self, other_state: Self) -> EulerCellDelta {
            EulerCellDelta {
                density: self.density - other_state.density,
                momentum: self.momentum - other_state.momentum,
                energy: self.energy - other_state.energy,
            }
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
        pub(crate) fn inner_flux(roe_avg: PrimitiveState, u_delta: EulerCellDelta, abs_eigenvalues: (f64, f64, f64)) -> EulerFlux {
            // left: EulerCell1d, right: EulerCell1d)
            // let roe_avg = Roe::average_state(left.to_primitive(), right.to_primitive());


//             dbg!(abs_eigenvalues);
//
//             dbg!(Roe::first_wave_strength(roe_avg, u_delta));
//             dbg!(Roe::first_eigvec(roe_avg));
//             dbg!(Roe::second_wave_strength(roe_avg, u_delta));
//             dbg!(Roe::second_eigvec(roe_avg));
//             dbg!(Roe::third_wave_strength(roe_avg, u_delta));
//             dbg!(Roe::third_eigvec(roe_avg));
//
            // find wave strengths
            -0.5 * (abs_eigenvalues.0
                * Roe::first_wave_strength(roe_avg, u_delta)
                * Roe::first_eigvec(roe_avg)
                + abs_eigenvalues.1
                    * Roe::second_wave_strength(roe_avg, u_delta)
                    * Roe::second_eigvec(roe_avg)
                + abs_eigenvalues.2
                    * Roe::third_wave_strength(roe_avg, u_delta)
                    * Roe::third_eigvec(roe_avg))

            // unimplemented!();
        }

        // associated eigenvalues
        pub fn eigenvalues(roe_avg: PrimitiveState) -> (f64, f64, f64) {
            let e_vals = roe_avg.eigenvalues();
            (
                f64::abs(e_vals.0),
                f64::abs(e_vals.1),
                f64::abs(e_vals.2),
            )
        }
        // -- right eigenvectors (K)
        // we represent them as EulerFluxes

        // Toro's K1
        #[inline]
        pub(crate) fn first_eigvec(roe_avg: PrimitiveState) -> EulerFlux {
            EulerFlux {
                density_flux: 1.0,
                momentum_flux: roe_avg.velocity - roe_avg.sound_speed(),
                energy_flux: roe_avg.enthalpy() - roe_avg.velocity * roe_avg.sound_speed(),
            }
        }

        // Toro's K2
        #[inline]
        pub(crate) fn second_eigvec(roe_avg: PrimitiveState) -> EulerFlux {
            EulerFlux {
                density_flux: 1.0,
                momentum_flux: roe_avg.velocity,
                energy_flux: 0.5 * roe_avg.velocity * roe_avg.velocity,
            }
        }

        // Toro's K3
        #[inline]
        pub(crate) fn third_eigvec(roe_avg: PrimitiveState) -> EulerFlux {
            EulerFlux {
                density_flux: 1.0,
                momentum_flux: roe_avg.velocity + roe_avg.sound_speed(),
                energy_flux: roe_avg.enthalpy() + roe_avg.velocity * roe_avg.sound_speed(),
            }
        }

        // -- wave strength coefficients (alpha)

        #[inline]
        pub(crate) fn first_wave_strength(roe_avg: PrimitiveState, u_delta: EulerCellDelta) -> f64 {
            1.0 / (2.0 * roe_avg.sound_speed())
                * (u_delta.density * (roe_avg.velocity + roe_avg.sound_speed())
                    - u_delta.momentum
                    - roe_avg.sound_speed() * Self::second_wave_strength(roe_avg, u_delta))
        }

        #[inline]
        pub(crate) fn second_wave_strength(
            roe_avg: PrimitiveState,
            u_delta: EulerCellDelta,
        ) -> f64 {
            (roe_avg.gamma - 1.0) / (roe_avg.sound_speed() * roe_avg.sound_speed())
                * (u_delta.density * (roe_avg.enthalpy() - roe_avg.velocity * roe_avg.velocity)
                    + roe_avg.velocity * u_delta.momentum
                    - u_delta.energy)
        }

        #[inline]
        pub(crate) fn third_wave_strength(roe_avg: PrimitiveState, u_delta: EulerCellDelta) -> f64 {
            u_delta.density
                - (Self::first_wave_strength(roe_avg, u_delta)
                    + Self::second_wave_strength(roe_avg, u_delta))
        }
    }

    #[derive(Debug, Clone, Copy)]
    pub(crate) struct EulerCellDelta {
        pub density: f64,
        pub momentum: f64,
        pub energy: f64,
    }

    // Roe with entropy fix takes advantage of the Roe functions
    #[derive(Debug, Clone, Copy)]
    pub struct RoeEntropyFix;

    impl FluxFunction for RoeEntropyFix {
        fn calculate_flux(
            &self,
            left: EulerCell1d,
            right: EulerCell1d,
            _time_step: f64,
        ) -> EulerFlux {

           let roe_avg = Roe::average_state(left.to_primitive(), right.to_primitive());

            0.5 * (left.to_flux()
                    + right.to_flux())
                + Roe::inner_flux(
                    roe_avg,
                    right.state_delta(left),
                    RoeEntropyFix::eigenvalues(left.to_primitive(), right.to_primitive()),
                )
        }
    }

    impl RoeEntropyFix {
        pub fn eigenvalues(left_state: PrimitiveState, right_state: PrimitiveState) -> (f64, f64, f64) {
            let left_eigvals = left_state.eigenvalues();
            let right_eigvals = right_state.eigenvalues();

            let abs_avg_eigvals = Roe::eigenvalues(Roe::average_state(left_state, right_state));

            (
                // could be a lot nicer with mapping over array, slice etc
                Self::harten_entropy_fix(left_eigvals.0, right_eigvals.0, abs_avg_eigvals.0),
                Self::harten_entropy_fix(left_eigvals.1, right_eigvals.1, abs_avg_eigvals.1),
                Self::harten_entropy_fix(left_eigvals.2, right_eigvals.2, abs_avg_eigvals.2),
            )
        }

        pub (crate) fn harten_entropy_fix(l_val: f64, r_val: f64, abs_avg_val: f64) -> f64 {
            // shock -- no problem, return roe_avg
            let scaling_factor = Self::eig_scaling_factor(l_val, r_val);

            // possible issue with div 0 -- check that scaling_factor isn't too small
            if l_val > r_val || scaling_factor < 1e-6 {
                return abs_avg_val
            }

            if abs_avg_val > scaling_factor / 2.0 {
                abs_avg_val
            } else {
                abs_avg_val * abs_avg_val / scaling_factor + scaling_factor / 4.0
            }

        }

        // little delta scaling factor
        #[inline]
        pub (crate) fn eig_scaling_factor(l_val: f64, r_val: f64) -> f64 {
            f64::max(0.0, 4.0 * (r_val - l_val))
        }
    }


    #[cfg(test)]
    mod tests {
        use super::*;
        use approx::*;

        #[test]
        fn roe_average() {
            let left_state = PrimitiveState {
                density: 2.281,
                velocity: 164.83,
                pressure: 201.17e3,
                gamma: 1.4,
            };

            assert_relative_eq!(
                Roe::average_density(left_state, left_state.clone()),
                left_state.density,
                epsilon = 1e-6
            );
            assert_relative_eq!(
                Roe::average_velocity(left_state, left_state.clone()),
                left_state.velocity,
                epsilon = 1e-6
            );
            assert_relative_eq!(
                Roe::average_pressure(left_state, left_state.clone()),
                left_state.pressure,
                epsilon = 1e-6
            );
            assert_relative_eq!(
                Roe::average_enthalpy(left_state, left_state.clone()),
                left_state.enthalpy(),
                epsilon = 1e-6
            );
        }

        #[test]
        fn inner_flux_direction() {
            let left_cell = EulerCell1d::new(
                PrimitiveState {
                    density: 2.281,
                    velocity: 164.83,
                    pressure: 201.17e3,
                    gamma: 1.4,
                },
                (0.0, 1.0)
            );

            let right_cell = EulerCell1d::new(
                PrimitiveState {
                    density: 1.408,
                    velocity: 0.0,
                    pressure: 101.1e3,
                    gamma: 1.4,
                },
                (1.0, 2.0));

            let roe_avg = Roe::average_state(left_cell.to_primitive(), right_cell.to_primitive());
            let roe_avg_flux = Roe::inner_flux(roe_avg, right_cell.state_delta(left_cell), Roe::eigenvalues(roe_avg));

            // expect everything moving to the right
            assert!(roe_avg_flux.density_flux > 0.0);
            // not sure about this one
            assert!(roe_avg_flux.momentum_flux > 0.0);
            assert!(roe_avg_flux.energy_flux > 0.0);

            dbg!(roe_avg_flux);

            let right_cell = EulerCell1d::new(
                PrimitiveState {
                    density: 2.281,
                    velocity: -164.83,
                    pressure: 201.17e3,
                    gamma: 1.4,
                },
                (1.0, 2.0)
            );

            let left_cell = EulerCell1d::new(
                PrimitiveState {
                    density: 1.408,
                    velocity: 0.0,
                    pressure: 101.1e3,
                    gamma: 1.4,
                },
                (0.0, 1.0)
            );

            let roe_avg = Roe::average_state(left_cell.to_primitive(), right_cell.to_primitive());
            let roe_avg_flux = Roe::inner_flux(roe_avg, right_cell.state_delta(left_cell), Roe::eigenvalues(roe_avg));

            assert!(roe_avg_flux.density_flux < 0.0);
            // not sure about this one
            assert!(roe_avg_flux.momentum_flux > 0.0);
            assert!(roe_avg_flux.energy_flux < 0.0);

            dbg!(roe_avg_flux);

        }

        #[test]
        fn roe_flux_direction() {
            let left_cell = EulerCell1d::new(
            /*let left_state = */ PrimitiveState {
                density: 2.281,
                velocity: 164.83,
                pressure: 201.17e3,
                gamma: 1.4,
            },
            (0.0, 1.0)
            );

            let right_cell = EulerCell1d::new(
                PrimitiveState {
                    density: 1.408,
                    velocity: 0.0,
                    pressure: 101.1e3,
                    gamma: 1.4,
                },
                (1.0, 2.0));

            let roe_avg_flux = Roe.calculate_flux(left_cell, right_cell, 0.0);

            // expect everything moving to the right
            assert!(roe_avg_flux.density_flux > 0.0);
            assert!(roe_avg_flux.momentum_flux > 0.0);
            assert!(roe_avg_flux.energy_flux > 0.0);

            dbg!(roe_avg_flux);

            let right_cell = EulerCell1d::new(
                PrimitiveState {
                    density: 2.281,
                    velocity: -164.83,
                    pressure: 201.17e3,
                    gamma: 1.4,
                },
                (1.0, 2.0)
            );

            let left_cell = EulerCell1d::new(
                PrimitiveState {
                    density: 1.408,
                    velocity: 0.0,
                    pressure: 101.1e3,
                    gamma: 1.4,
                },
                (0.0, 1.0)
            );


            // test reverse
            let roe_avg_flux = Roe.calculate_flux(left_cell, right_cell, 0.0);
            assert!(roe_avg_flux.density_flux < 0.0);
            assert!(roe_avg_flux.momentum_flux > 0.0);
            assert!(roe_avg_flux.energy_flux < 0.0);

            dbg!(roe_avg_flux);
        }

        #[test]
        fn empty_roe_flux() {
            let left_state = PrimitiveState {
                density: 2.281,
                velocity: 164.83,
                pressure: 201.17e3,
                gamma: 1.4,
            };

            let zero_delta = EulerCellDelta {
                density: 0.0,
                momentum: 0.0,
                energy: 0.0,
            };

            let zero_flux = EulerFlux {
                density_flux: 0.0,
                momentum_flux: 0.0,
                energy_flux: 0.0,
            };

            let roe_avg_flux = Roe::inner_flux(left_state, zero_delta, Roe::eigenvalues(left_state));

            assert_relative_eq!(roe_avg_flux.density_flux, zero_flux.density_flux);
            assert_relative_eq!(roe_avg_flux.momentum_flux, zero_flux.momentum_flux);
            assert_relative_eq!(roe_avg_flux.energy_flux, zero_flux.energy_flux);
        }

        #[test]
        fn wave_strengths() {
            let left_state = PrimitiveState {
                density: 2.281,
                velocity: 164.83,
                pressure: 201.17e3,
                gamma: 1.4,
            };

            let cell_delta = EulerCellDelta {
                density: 0.0,
                momentum: 0.0,
                energy: 0.0,
            };

            let nonzero_delta = EulerCellDelta {
                density: 1.0,
                momentum: 1.0,
                energy: 1.0,
            };

            assert_relative_eq!(
                Roe::first_wave_strength(left_state, cell_delta),
                0.0,
                epsilon = 1e-6
            );

            assert_relative_ne!(
                Roe::first_wave_strength(left_state, nonzero_delta),
                0.0,
                epsilon = 1e-6
            );

            assert_relative_eq!(
                Roe::second_wave_strength(left_state, cell_delta),
                0.0,
                epsilon = 1e-6
            );

            assert_relative_ne!(
                Roe::second_wave_strength(left_state, nonzero_delta),
                0.0,
                epsilon = 1e-6
            );

            assert_relative_eq!(
                Roe::third_wave_strength(left_state, cell_delta),
                0.0,
                epsilon = 1e-6
            );

            assert_relative_ne!(
                Roe::third_wave_strength(left_state, nonzero_delta),
                0.0,
                epsilon = 1e-6
            );
        }
    }
}

pub mod second_order {
    use super::*;
}
