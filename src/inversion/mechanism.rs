//! Moment tensors from the base-element weights: the determinant routines,
//! the double-couple constraint `zmom` and the decomposition `dcp`. A tensor
//! is six weights `v`, the last isotropic, laid out as `mtrx` reads them.

use super::linalg::{self, Matrix, N};

pub type Weights = [f32; N];

/// The symmetric moment tensor, x north, y east, z down.
pub fn mtrx(v: &Weights) -> [[f32; 3]; 3] {
    let m11 = (v[1] - v[4]) + v[5];
    let m22 = -v[1] + v[5];
    let m33 = v[4] + v[5];
    [[m11, v[0], v[3]], [v[0], m22, v[2]], [v[3], v[2], m33]]
}

/// The tensor as nine reals in Fortran column order.
pub fn columns(v: &Weights) -> [f32; 9] {
    let m = mtrx(v);
    std::array::from_fn(|i| m[i % 3][i / 3])
}

/// The cofactors of the first row and the other three minors `det1` uses,
/// each as the reference groups it.
fn minors(m: &[[f32; 3]; 3]) -> [f32; 6] {
    [
        m[1][1] * m[2][2] - m[1][2] * m[1][2],
        m[0][1] * m[2][2] - m[1][2] * m[0][2],
        m[0][1] * m[1][2] - m[1][1] * m[0][2],
        m[0][0] * m[2][2] - m[0][2] * m[0][2],
        m[0][0] * m[1][2] - m[0][1] * m[0][2],
        m[0][0] * m[1][1] - m[0][1] * m[0][1],
    ]
}

/// The determinant, computed in single precision as the reference does.
pub fn det0(v: &Weights) -> f64 {
    let m = mtrx(v);
    let c = minors(&m);
    (m[0][0] * c[0] - m[0][1] * c[1] + m[0][2] * c[2]) as f64
}

/// The gradient of the determinant over the five deviatoric weights.
pub fn det00(v: &Weights) -> Weights {
    let [a11, a12, a13, a22, a23, a33] = minors(&mtrx(v));
    [
        -(2.0 * a12),
        a11 - a22,
        -(2.0 * a23),
        2.0 * a13,
        -a11 + a33,
        0.0,
    ]
}

/// The determinant with one factor of `v`'s tensor taken from `u`'s.
pub fn det1(u: &Weights, v: &Weights) -> f64 {
    let a = mtrx(u);
    let c = minors(&mtrx(v));
    (a[0][0] * c[0] - a[0][1] * c[1] * 2.0 + a[0][2] * c[2] * 2.0 + a[1][1] * c[3]
        - a[1][2] * c[4] * 2.0
        + a[2][2] * c[5]) as f64
}

/// `zmom`: the zero-determinant tensor nearest the unconstrained one. The
/// reference keeps the cubic's roots between calls, and a failed solve
/// reuses them, so they live here.
#[derive(Default)]
pub struct Zmom {
    xr: [f64; 3],
    xi: [f64; 3],
}

impl Zmom {
    /// Moves `r2` to the weights nearest `r1`, in the norm of `rww`, whose
    /// tensor has zero determinant; `arw` is the inverse of `rww`. Returns
    /// the solver's flag, 0 when the last cubic was solved, or `None` for a
    /// zero `r1`, where `r2` is left alone.
    pub fn solve(
        &mut self,
        arw: &Matrix,
        rww: &Matrix,
        r1: &Weights,
        r2: &mut Weights,
    ) -> Option<i32> {
        let mut dob0 = 0f32;
        let mut an = 0f32;
        for n1 in 0..N {
            for n2 in 0..N {
                an += rww[n1][n2] * r1[n1] * r1[n2];
            }
        }
        if an == 0.0 {
            return None;
        }
        let mut da = [0f64; 4];
        da[3] = det0(r1);
        let mut uj = [0f32; N];
        let mut rr = [0f32; N];
        let mut il = 0;
        for i in 1..=30 {
            let u = det00(r2);
            let mut ui = [0f32; N];
            for n in 0..N {
                for n1 in 0..N {
                    ui[n] -= arw[n][n1] * u[n1];
                }
                if i > 1 {
                    ui[n] = (ui[n] + uj[n]) / 2.0;
                }
                uj[n] = ui[n];
            }
            da[2] = det1(&ui, r1);
            da[1] = det1(r1, &ui);
            da[0] = det0(&ui);
            il = if linalg::cubic(da, &mut self.xr, &mut self.xi) {
                0
            } else {
                30000
            };
            let mut al = self.xr[0];
            for k in 1..3 {
                if self.xi[k] == 0.0 && !(self.xr[k].abs() >= al.abs()) {
                    al = self.xr[k];
                }
            }
            for n in 0..N {
                r2[n] = (r1[n] as f64 + al * ui[n] as f64) as f32;
            }
            let mut dob = 0f32;
            for n1 in 0..N {
                for n2 in 0..N {
                    dob += rww[n1][n2] * (r2[n1] - r1[n1]) * (r2[n2] - r1[n2]);
                }
            }
            if i > 3 && dob >= dob0 {
                continue;
            }
            rr = *r2;
            if (dob - dob0).abs() / an < 1e-4 {
                break;
            }
            dob0 = dob;
        }
        *r2 = rr;
        Some(il)
    }
}

/// A tensor read as a double couple.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DoubleCouple {
    /// Half the spread of the largest and smallest eigenvalues.
    pub moment: f32,
    /// `moment` less the largest eigenvalue.
    pub non_dc: f32,
    /// The fault plane, in degrees.
    pub strike: f32,
    pub dip: f32,
    pub slip: f32,
    /// The tension axis.
    pub axis: [f32; 3],
}

/// `dcp`. `None` when the eigenvalues do not converge, where the reference
/// pauses.
pub fn dcp(v: &Weights) -> Option<DoubleCouple> {
    const PI: f32 = 3.1415926;
    use std::f32::consts::SQRT_2 as ROOT2;
    let rad = PI / 180.0;
    let (mut eg, mut ev) = linalg::jacobi(mtrx(v))?;
    linalg::eigsrt(&mut eg, &mut ev);
    let (mut e1, mut e3, mut i1, mut i3) = (eg[0], eg[0], 0, 0);
    for ii in 1..3 {
        if !(e1 >= eg[ii]) {
            e1 = eg[ii];
            i1 = ii;
        }
        if !(e3 < eg[ii]) {
            e3 = eg[ii];
            i3 = ii;
        }
    }
    let moment = (e1 - e3) / 2.0;
    let non_dc = moment - e1;
    let axis = [ev[0][i1], ev[1][i1], ev[2][i1]];
    let slip_vector: [f32; 3] = std::array::from_fn(|i| (ev[i][i3] - ev[i][i1]) / ROOT2);
    let normal: [f32; 3] = std::array::from_fn(|i| (ev[i][i1] + ev[i][i3]) / ROOT2);
    let mut dip = normal[2].acos();
    let mut strike = normal[1].atan2(normal[0]) + PI / 2.0;
    let (cf, sf) = (strike.cos(), strike.sin());
    let c = slip_vector[0] * cf + slip_vector[1] * sf;
    let s = if dip != 0.0 {
        slip_vector[2] / dip.sin()
    } else {
        slip_vector[0] * sf - slip_vector[1] * cf
    };
    strike /= rad;
    dip /= rad;
    let mut slip = -s.atan2(c) / rad;
    if !(dip < 90.0) {
        slip = -slip;
        dip = 180.0 - dip;
        strike += 180.0;
        if strike > 360.0 {
            strike -= 360.0;
        }
    }
    Some(DoubleCouple {
        moment,
        non_dc,
        strike,
        dip,
        slip,
        axis,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dcp_reads_a_vertical_strike_slip() {
        // Mxy = 1: a vertical strike-slip of moment 1.
        let mechanism = dcp(&[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]).unwrap();
        assert!((mechanism.moment - 1.0).abs() < 1e-6);
        assert!(mechanism.non_dc.abs() < 1e-6);
        assert!((mechanism.dip - 90.0).abs() < 1e-3);
        assert!(mechanism.slip.abs() % 180.0 < 1e-3);
    }

    #[test]
    fn zmom_gives_a_tensor_with_zero_determinant() {
        let mut identity = [[0f32; N]; N];
        for (n, row) in identity.iter_mut().enumerate().take(5) {
            row[n] = 1.0;
        }
        let r1 = [1.0, 0.3, 0.2, 0.1, 0.4, 0.0];
        let mut r2 = r1;
        let il = Zmom::default().solve(&identity, &identity, &r1, &mut r2);
        assert_eq!(il, Some(0));
        assert!(det0(&r2).abs() < 1e-4 * det0(&r1).abs());
    }
}
