//! The Numerical Recipes routines `inversion` links (`ludcmp`, `lubksb`,
//! `jacobi`, `eigsrt`) and the Nagoya University cubic and quadratic solvers.
//! Each keeps the reference's operation order, so results match it bit for
//! bit. Matrices are indexed `[row][column]`, as `a(row, column)`.

/// The largest matrix: six base elements.
pub const N: usize = 6;

pub type Matrix = [[f32; N]; N];

/// Crout LU decomposition with partial pivoting of the leading `n` by `n`
/// block of `a`, in place. Returns the row permutation, or `None` for a
/// matrix with a zero row, where the reference pauses.
pub fn ludcmp(a: &mut Matrix, n: usize) -> Option<[usize; N]> {
    let mut vv = [0f32; N];
    for i in 0..n {
        let mut aamax = 0f32;
        for j in 0..n {
            if a[i][j].abs() > aamax {
                aamax = a[i][j].abs();
            }
        }
        if aamax == 0.0 {
            return None;
        }
        vv[i] = 1.0 / aamax;
    }
    let mut indx = [0; N];
    let mut imax = 0;
    for j in 0..n {
        for i in 0..j {
            let mut sum = a[i][j];
            for k in 0..i {
                sum -= a[i][k] * a[k][j];
            }
            a[i][j] = sum;
        }
        let mut aamax = 0f32;
        for i in j..n {
            let mut sum = a[i][j];
            for k in 0..j {
                sum -= a[i][k] * a[k][j];
            }
            a[i][j] = sum;
            let dum = vv[i] * sum.abs();
            if dum >= aamax {
                imax = i;
                aamax = dum;
            }
        }
        if j != imax {
            for k in 0..n {
                let dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if a[j][j] == 0.0 {
            a[j][j] = 1.0e-20;
        }
        if j != n - 1 {
            let dum = 1.0 / a[j][j];
            for row in a.iter_mut().take(n).skip(j + 1) {
                row[j] *= dum;
            }
        }
    }
    Some(indx)
}

/// Solves `A x = b` in place, with `a` and `indx` from [`ludcmp`].
pub fn lubksb(a: &Matrix, n: usize, indx: &[usize; N], b: &mut [f32; N]) {
    let mut first = None;
    for i in 0..n {
        let ll = indx[i];
        let mut sum = b[ll];
        b[ll] = b[i];
        if let Some(ii) = first {
            for j in ii..i {
                sum -= a[i][j] * b[j];
            }
        } else if sum != 0.0 {
            first = Some(i);
        }
        b[i] = sum;
    }
    for i in (0..n).rev() {
        let mut sum = b[i];
        for j in i + 1..n {
            sum -= a[i][j] * b[j];
        }
        b[i] = sum / a[i][i];
    }
}

/// Eigenvalues and eigenvectors (as columns) of a symmetric 3 by 3 matrix,
/// by Jacobi rotations on its upper triangle. `None` after 50 sweeps, where
/// the reference pauses.
pub fn jacobi(mut a: [[f32; 3]; 3]) -> Option<([f32; 3], [[f32; 3]; 3])> {
    const SIZE: usize = 3;
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut b = [a[0][0], a[1][1], a[2][2]];
    let mut d = b;
    let mut z = [0f32; SIZE];
    for sweep in 1..=50 {
        let mut sm = 0f32;
        for p in 0..SIZE - 1 {
            for q in p + 1..SIZE {
                sm += a[p][q].abs();
            }
        }
        if sm == 0.0 {
            return Some((d, v));
        }
        let tresh = if sweep < 4 { 0.2 * sm / 9.0 } else { 0.0 };
        for p in 0..SIZE - 1 {
            for q in p + 1..SIZE {
                let g = 100.0 * a[p][q].abs();
                if sweep > 4 && d[p].abs() + g == d[p].abs() && d[q].abs() + g == d[q].abs() {
                    a[p][q] = 0.0;
                } else if a[p][q].abs() > tresh {
                    let h = d[q] - d[p];
                    let t = if h.abs() + g == h.abs() {
                        a[p][q] / h
                    } else {
                        let theta = 0.5 * h / a[p][q];
                        let t = 1.0 / (theta.abs() + (1.0 + theta * theta).sqrt());
                        if theta < 0.0 { -t } else { t }
                    };
                    let c = 1.0 / (1.0 + t * t).sqrt();
                    let s = t * c;
                    let tau = s / (1.0 + c);
                    let h = t * a[p][q];
                    z[p] -= h;
                    z[q] += h;
                    d[p] -= h;
                    d[q] += h;
                    a[p][q] = 0.0;
                    let rotate = |g: f32, h: f32| (g - s * (h + g * tau), h + s * (g - h * tau));
                    for j in 0..p {
                        (a[j][p], a[j][q]) = rotate(a[j][p], a[j][q]);
                    }
                    for j in p + 1..q {
                        (a[p][j], a[j][q]) = rotate(a[p][j], a[j][q]);
                    }
                    for j in q + 1..SIZE {
                        (a[p][j], a[q][j]) = rotate(a[p][j], a[q][j]);
                    }
                    for row in &mut v {
                        (row[p], row[q]) = rotate(row[p], row[q]);
                    }
                }
            }
        }
        for p in 0..SIZE {
            b[p] += z[p];
            d[p] = b[p];
            z[p] = 0.0;
        }
    }
    None
}

/// Sorts eigenvalues into descending order by straight insertion, moving
/// the eigenvector columns with them.
pub fn eigsrt(d: &mut [f32; 3], v: &mut [[f32; 3]; 3]) {
    for i in 0..2 {
        let mut k = i;
        let mut p = d[i];
        for (j, &value) in d.iter().enumerate().skip(i + 1) {
            if value >= p {
                k = j;
                p = value;
            }
        }
        if k != i {
            d[k] = d[i];
            d[i] = p;
            for row in v.iter_mut() {
                row.swap(i, k);
            }
        }
    }
}

/// The roots `x + iy` of `a[0] t³ + a[1] t² + a[2] t + a[3]`, the real root
/// first. Returns `false` when `a[0]` is zero, leaving `x` and `y` as they
/// were. A zero leading term of the deflated quadratic leaves its two roots
/// as they were too.
pub fn cubic(a: [f64; 4], x: &mut [f64; 3], y: &mut [f64; 3]) -> bool {
    let mut w = a;
    if w[0] == 0.0 {
        return false;
    }
    y[0] = 0.0;
    if w[3] == 0.0 {
        x[0] = 0.0;
    } else {
        let c = w[0] * 3.0;
        let b = w[2] / c;
        let c = w[1] / c;
        let p = c * c - b;
        let mut q = (0.5 * b - p) * c - w[3] * 0.5 / w[0];
        if q != 0.0 {
            let d = q * q - p * p * p;
            q = if d < 0.0 {
                ((-d).sqrt().atan2(q) / 3.0).cos() * p.sqrt() * 2.0
            } else if d == 0.0 {
                (p.sqrt() * 2.0).copysign(q)
            } else {
                let d1 = d.sqrt().copysign(q) + q;
                let u = d1.abs().powf(1.0 / 3.0).copysign(d1);
                let v = p / u;
                if p < 0.0 {
                    (q + q) / ((u - v) * u + v * v)
                } else {
                    u + v
                }
            };
        }
        x[0] = q - c;
        if x[0].abs() <= w[3].abs() / 3.0 {
            w[1] = x[0] * w[0] + w[1];
            w[2] = x[0] * w[1] + w[2];
        } else {
            let p = w[1];
            w[1] = w[3] / x[0] + w[2];
            w[0] = w[1] / x[0] + p;
            w[2] = w[3];
        }
    }
    let (_, xs) = x.split_at_mut(1);
    let (_, ys) = y.split_at_mut(1);
    quadratic([w[0], w[1], w[2]], xs, ys);
    true
}

/// The roots `x + iy` of `a[0] t² + a[1] t + a[2]`. Leaves `x` and `y` as
/// they were when `a[0]` is zero.
fn quadratic(a: [f64; 3], x: &mut [f64], y: &mut [f64]) {
    if a[0] == 0.0 {
        return;
    }
    if a[2] == 0.0 {
        x[0] = -a[1] / a[0];
        x[1] = 0.0;
        y[0] = 0.0;
        y[1] = 0.0;
        return;
    }
    let p = -a[1] / a[0] * 0.5;
    let q = a[2] / a[0];
    let d = p * p - q;
    if d < 0.0 {
        y[0] = (-d).sqrt();
        y[1] = -y[0];
        x[0] = p;
        x[1] = p;
    } else {
        x[0] = d.sqrt().copysign(p) + p;
        x[1] = q / x[0];
        y[0] = 0.0;
        y[1] = 0.0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ludcmp_and_lubksb_invert_a_matrix() {
        let mut a = [[0f32; N]; N];
        a[0][..3].copy_from_slice(&[2.0, 1.0, 0.0]);
        a[1][..3].copy_from_slice(&[1.0, 3.0, 1.0]);
        a[2][..3].copy_from_slice(&[0.0, 1.0, 4.0]);
        let original = a;
        let indx = ludcmp(&mut a, 3).unwrap();
        let mut b = [1.0, 2.0, 3.0, 0.0, 0.0, 0.0];
        lubksb(&a, 3, &indx, &mut b);
        for row in &original[..3] {
            let product: f32 = row[..3].iter().zip(&b).map(|(m, x)| m * x).sum();
            assert!(
                (product - [1.0, 2.0, 3.0][original.iter().position(|r| r == row).unwrap()]).abs()
                    < 1e-5
            );
        }
    }

    #[test]
    fn ludcmp_refuses_a_zero_row() {
        let mut a = [[0f32; N]; N];
        a[0][0] = 1.0;
        assert_eq!(ludcmp(&mut a, 2), None);
    }

    #[test]
    fn jacobi_and_eigsrt_give_descending_eigenvalues() {
        let (mut d, mut v) = jacobi([[2.0, 1.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 5.0]]).unwrap();
        eigsrt(&mut d, &mut v);
        assert!(
            (d[0] - 5.0).abs() < 1e-6 && (d[1] - 3.0).abs() < 1e-6 && (d[2] - 1.0).abs() < 1e-6
        );
        assert!((v[2][0].abs() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn cubic_finds_three_real_roots() {
        // (t - 1)(t - 2)(t + 3) = t³ - 7t + 6
        let (mut x, mut y) = ([0.0; 3], [0.0; 3]);
        assert!(cubic([1.0, 0.0, -7.0, 6.0], &mut x, &mut y));
        let mut roots = x;
        roots.sort_by(f64::total_cmp);
        for (root, expected) in roots.iter().zip([-3.0, 1.0, 2.0]) {
            assert!((root - expected).abs() < 1e-12);
        }
        assert_eq!(y, [0.0; 3]);
    }

    #[test]
    fn cubic_leaves_the_roots_when_it_is_not_cubic() {
        let (mut x, mut y) = ([7.0; 3], [8.0; 3]);
        assert!(!cubic([0.0, 1.0, 2.0, 3.0], &mut x, &mut y));
        assert_eq!((x, y), ([7.0; 3], [8.0; 3]));
    }
}
