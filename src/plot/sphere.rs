//! The focal sphere on a lower-hemisphere equal-area net, as the reference's
//! `stnplo` and `mcplot0` draw it: x north, y east, z down; on the page, east
//! is right and north up.

use crate::formats::solution::Mechanism;

pub type Vector = [f64; 3];

/// The tensor from nine reals in Fortran column order.
pub fn tensor(columns: &[f32; 9]) -> [[f64; 3]; 3] {
    std::array::from_fn(|i| std::array::from_fn(|j| columns[i + 3 * j] as f64))
}

/// The P-wave amplitude `dᵀMd` along the ray `d`, `radpmt`'s `pw`.
pub fn p_amplitude(m: &[[f64; 3]; 3], d: Vector) -> f64 {
    (0..3)
        .flat_map(|i| (0..3).map(move |j| (i, j)))
        .map(|(i, j)| m[i][j] * d[i] * d[j])
        .sum()
}

/// The ray leaving at `azimuth` from north and `takeoff` from down, both in
/// radians.
pub fn ray(azimuth: f64, takeoff: f64) -> Vector {
    let (st, ct) = takeoff.sin_cos();
    let (sa, ca) = azimuth.sin_cos();
    [st * ca, st * sa, ct]
}

/// The point on the net of radius 1 for a ray, folded onto the lower
/// hemisphere, as (right, up).
pub fn project(d: Vector) -> (f64, f64) {
    let d = if d[2] < 0.0 { d.map(|c| -c) } else { d };
    let takeoff = d[2].clamp(-1.0, 1.0).acos();
    let r = std::f64::consts::SQRT_2 * (takeoff / 2.0).sin();
    let azimuth = d[1].atan2(d[0]);
    (r * azimuth.sin(), r * azimuth.cos())
}

/// The ray under a point of the net of radius 1, the inverse of `project`.
pub fn unproject(right: f64, up: f64) -> Vector {
    let r = right.hypot(up).min(1.0);
    let takeoff = 2.0 * (r / std::f64::consts::SQRT_2).asin();
    ray(right.atan2(up), takeoff)
}

/// The fault normal and the slip vector, Aki and Richards' (4.88).
pub fn normal_and_slip(m: Mechanism) -> (Vector, Vector) {
    let (ss, cs) = (m.strike as f64).to_radians().sin_cos();
    let (sd, cd) = (m.dip as f64).to_radians().sin_cos();
    let (sl, cl) = (m.slip as f64).to_radians().sin_cos();
    let normal = [-sd * ss, sd * cs, -cd];
    let slip = [cl * cs + cd * sl * ss, cl * ss - cd * sl * cs, -sl * sd];
    (normal, slip)
}

/// The lower half of the great circle perpendicular to `pole`, in `n` + 1
/// points on the net of radius 1.
pub fn great_circle(pole: Vector, n: usize) -> Vec<(f64, f64)> {
    let horizontal = pole[0].hypot(pole[1]);
    let e1 = if horizontal < 1e-9 {
        [1.0, 0.0, 0.0]
    } else {
        [-pole[1] / horizontal, pole[0] / horizontal, 0.0]
    };
    let mut e2 = cross(pole, e1);
    if e2[2] < 0.0 {
        e2 = e2.map(|c| -c);
    }
    let norm = e2.iter().map(|c| c * c).sum::<f64>().sqrt();
    (0..=n)
        .map(|k| {
            let t = std::f64::consts::PI * k as f64 / n as f64;
            let (s, c) = t.sin_cos();
            project(std::array::from_fn(|i| c * e1[i] + s * e2[i] / norm))
        })
        .collect()
}

/// The compressional part of the net of radius 1 as horizontal runs in `n`
/// rows: (up, left, right) for each run.
pub fn compressions(m: &[[f64; 3]; 3], n: usize) -> Vec<(f64, f64, f64)> {
    let step = 2.0 / n as f64;
    let mut runs = Vec::new();
    for row in 0..n {
        let up = 1.0 - (row as f64 + 0.5) * step;
        let half = (1.0 - up * up).sqrt();
        let columns = ((2.0 * half / step).ceil() as usize).max(1);
        let width = 2.0 * half / columns as f64;
        let mut start = None;
        for column in 0..=columns {
            let inside = column < columns && {
                let right = -half + (column as f64 + 0.5) * width;
                p_amplitude(m, unproject(right, up)) > 0.0
            };
            let left = -half + column as f64 * width;
            match (inside, start) {
                (true, None) => start = Some(left),
                (false, Some(s)) => {
                    runs.push((up, s, left));
                    start = None;
                }
                _ => {}
            }
        }
    }
    runs
}

pub fn cross(a: Vector, b: Vector) -> Vector {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
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

    /// The double couple `n uᵀ + u nᵀ`.
    fn couple(m: Mechanism) -> [[f64; 3]; 3] {
        let (n, u) = normal_and_slip(m);
        std::array::from_fn(|i| std::array::from_fn(|j| n[i] * u[j] + u[i] * n[j]))
    }

    #[test]
    fn a_thrust_compresses_the_centre_and_dilates_the_rim() {
        let m = couple(mechanism(0.0, 45.0, 90.0));
        assert!(p_amplitude(&m, unproject(0.0, 0.0)) > 0.0);
        assert!(p_amplitude(&m, unproject(0.99, 0.0)) < 0.0);
    }

    #[test]
    fn projection_round_trips() {
        let d = ray(1.0, 0.7);
        let (right, up) = project(d);
        let back = unproject(right, up);
        assert!((0..3).all(|i| (back[i] - d[i]).abs() < 1e-12));
    }

    #[test]
    fn both_nodal_planes_lie_on_zero_amplitude() {
        let m = mechanism(310.0, 20.0, 75.0);
        let (normal, slip) = normal_and_slip(m);
        let tensor = couple(m);
        for pole in [normal, slip] {
            for (right, up) in great_circle(pole, 50) {
                assert!(p_amplitude(&tensor, unproject(right, up)).abs() < 1e-9);
            }
        }
    }
}
