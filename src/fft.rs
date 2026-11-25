//! Single-precision complex arithmetic and the FFT of `sub.cfft.f`.
//!
//! The operations follow GCC's lowering of Fortran `COMPLEX` arithmetic at
//! `-O`: a full product, Smith's division, and the shorter forms GCC emits
//! when one operand is known to be real or imaginary. Keeping the same form
//! gives the reference's bits.

use std::ops::{Add, Div, Mul, Neg, Sub};

#[expect(
    clippy::approx_constant,
    reason = "the reference's value, one ulp above f32::consts::PI"
)]
pub const PI: f32 = 3.141593;

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Complex {
    pub re: f32,
    pub im: f32,
}

impl Complex {
    pub const fn new(re: f32, im: f32) -> Self {
        Self { re, im }
    }

    pub fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }

    /// The product with a real.
    pub fn scale(self, r: f32) -> Self {
        Self::new(r * self.re, r * self.im)
    }

    /// The product with the imaginary number `i c`.
    pub fn mul_imag(self, c: f32) -> Self {
        Self::new(-(self.im * c), self.re * c)
    }

    /// `exp(z)` as glibc's `cexpf` computes it for a finite argument whose
    /// real part is below 88, the only arguments the programs pass.
    pub fn exp(self) -> Self {
        let (sin, cos) = if self.im.abs() > f32::MIN_POSITIVE {
            self.im.sin_cos()
        } else {
            (self.im, 1.0)
        };
        let e = self.re.exp();
        Self::new(e * cos, e * sin)
    }

    /// `|z|`: `cabsf` is `hypotf`.
    pub fn abs(self) -> f32 {
        self.re.hypot(self.im)
    }

    /// `z ** n` by libgfortran's `pow_c4_i4`: binary powering from 1.
    pub fn powi(self, n: i32) -> Self {
        let mut pow = Self::new(1.0, 0.0);
        if n == 0 {
            return pow;
        }
        let mut x = if n < 0 { pow / self } else { self };
        let mut u = n.unsigned_abs();
        loop {
            if u & 1 == 1 {
                pow = pow * x;
            }
            u >>= 1;
            if u == 0 {
                return pow;
            }
            x = x * x;
        }
    }

    /// `log(i a)` for `a > 0`, by the branches glibc's `clogf` takes when the
    /// real part is zero.
    pub fn ln_imag(a: f32) -> Self {
        debug_assert!(a > 0.0 && a.is_finite());
        let re = if a == 1.0 {
            0.0
        } else if (0.5..2.0).contains(&a) {
            ((a - 1.0) * (a + 1.0)).ln_1p() / 2.0
        } else {
            a.ln()
        };
        Self::new(re, a.atan2(0.0))
    }
}

impl From<f32> for Complex {
    /// Fortran's conversion of a real operand in a complex expression.
    fn from(re: f32) -> Self {
        Self::new(re, 0.0)
    }
}

impl Neg for Complex {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.re, -self.im)
    }
}

/// `i w - b`, as GCC computes it with a zero real part on the left.
pub fn imag_minus(w: f32, b: Complex) -> Complex {
    Complex::new(-b.re, w - b.im)
}

impl Add for Complex {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self::new(self.re + b.re, self.im + b.im)
    }
}

impl Sub for Complex {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        Self::new(self.re - b.re, self.im - b.im)
    }
}

impl Mul for Complex {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self::new(
            self.re * b.re - self.im * b.im,
            self.re * b.im + self.im * b.re,
        )
    }
}

impl Div for Complex {
    type Output = Self;
    /// Smith's method, branching on the larger part of the divisor.
    fn div(self, b: Self) -> Self {
        let (ar, ai) = (self.re, self.im);
        if b.re.abs() < b.im.abs() {
            let ratio = b.re / b.im;
            let div = b.re * ratio + b.im;
            Self::new((ar * ratio + ai) / div, (ai * ratio - ar) / div)
        } else {
            let ratio = b.im / b.re;
            let div = b.im * ratio + b.re;
            Self::new((ai * ratio + ar) / div, (ai - ar * ratio) / div)
        }
    }
}

/// The unnormalised transform of `x` in place: `exp(-i...)` forward,
/// `exp(+i...)` inverse. The length is a power of two.
pub fn cfft(x: &mut [Complex], inverse: bool) {
    let n = x.len();
    let pi = if inverse { PI } else { -PI };
    let mut l = n;
    while l > 1 {
        let ni = l / 2;
        let angle = pi / ni as f32;
        let z = Complex::new(angle.cos(), angle.sin());
        let mut zk = Complex::new(1.0, 0.0);
        for k in 0..ni {
            for m in 0..n / l {
                let j = k + 2 * m * ni;
                let xj = x[j];
                x[j] = xj + x[j + ni];
                x[j + ni] = (xj - x[j + ni]) * zk;
            }
            zk = zk * z;
        }
        l = ni;
    }
    for i in 1..=n {
        let (mut j, mut m, mut l) = (i - 1, n / 2, 1);
        loop {
            l += m * (j % 2);
            m /= 2;
            j /= 2;
            if j < 1 {
                break;
            }
        }
        if i <= l {
            x.swap(i - 1, l - 1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_impulse_transforms_to_a_flat_spectrum_and_back() {
        let mut x = vec![Complex::default(); 8];
        x[1] = Complex::new(1.0, 0.0);
        cfft(&mut x, false);
        for (k, v) in x.iter().enumerate() {
            let angle = -std::f32::consts::PI * k as f32 / 4.0;
            assert!((v.re - angle.cos()).abs() < 1e-6 && (v.im - angle.sin()).abs() < 1e-6);
        }
        cfft(&mut x, true);
        assert!((x[1].re - 8.0).abs() < 1e-5);
        assert!(
            x.iter()
                .enumerate()
                .all(|(k, v)| k == 1 || v.re.abs() < 1e-5)
        );
    }
}
