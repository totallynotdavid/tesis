//! How far an inverted mechanism is from a catalogue one.

use crate::formats::solution::Mechanism;
use crate::plot::sphere::{Vector, cross, normal_and_slip};

/// Kagan's angle, in degrees: the smallest rotation that takes one double
/// couple onto the other. It is 0 for the same mechanism given by either
/// nodal plane, and at most 120.
pub fn kagan(a: Mechanism, b: Mechanism) -> f64 {
    let (ta, pa, ba) = axes(a);
    let (tb, pb, bb) = axes(b);
    let dot = |u: Vector, v: Vector| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    // A double couple is unchanged by a half turn about any of its axes,
    // which flips the sign of the other two.
    [
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ]
    .iter()
    .map(|s: &[f64; 3]| s[0] * dot(ta, tb) + s[1] * dot(pa, pb) + s[2] * dot(ba, bb))
    .map(|trace| ((trace - 1.0) / 2.0).clamp(-1.0, 1.0).acos().to_degrees())
    .fold(f64::INFINITY, f64::min)
}

/// The T, P and B axes.
fn axes(m: Mechanism) -> (Vector, Vector, Vector) {
    let (n, u) = normal_and_slip(m);
    let half = std::f64::consts::FRAC_1_SQRT_2;
    let t = std::array::from_fn(|i| (n[i] + u[i]) * half);
    let p = std::array::from_fn(|i| (n[i] - u[i]) * half);
    (t, p, cross(t, p))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mechanism(strike: f32, dip: f32, slip: f32) -> Mechanism {
        Mechanism {
            moment: 1.0,
            second: 0.0,
            strike,
            dip,
            slip,
        }
    }

    #[test]
    fn either_nodal_plane_names_the_same_double_couple() {
        // USGS's two planes for the 2013 Yauca Mwc solution.
        let angle = kagan(mechanism(131.0, 57.0, 90.0), mechanism(310.0, 33.0, 89.08));
        assert!(angle < 1.5, "{angle}");
    }

    #[test]
    fn a_vertical_strike_slip_turned_about_the_vertical_rotates_by_that_angle() {
        let angle = kagan(mechanism(10.0, 90.0, 0.0), mechanism(40.0, 90.0, 0.0));
        assert!((angle - 30.0).abs() < 1e-6, "{angle}");
    }

    #[test]
    fn swapping_the_t_and_p_axes_is_a_quarter_turn() {
        let angle = kagan(mechanism(130.0, 45.0, 90.0), mechanism(130.0, 45.0, -90.0));
        assert!((angle - 90.0).abs() < 1e-6, "{angle}");
    }
}
