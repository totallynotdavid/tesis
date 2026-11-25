//! Frequency responses on the half spectrum `0..nt/2`: the Q filter and the
//! two instrument types.

use super::c;
use crate::fft::{Complex, PI, imag_minus};

/// The attenuation filter `exp(i 2f t* log(i f/s2))` of `sub.qfm.f`, with
/// travel time over Q `tq`.
pub fn qfm(nt: usize, tq: f32, df: f32) -> Vec<Complex> {
    let s2 = 3.0 * ((2.0 * PI) / tq);
    let mut z = vec![c(1.0)];
    for i in 1..nt / 2 {
        let f = df * i as f32;
        z.push((Complex::new(0.0, (2.0 * f) * tq) * Complex::ln_imag(f / s2)).exp());
    }
    z
}

/// A response in poles and zeros, times `(iω)^power`, for a delta input.
pub fn instg(
    nt: usize,
    dw: f32,
    poles: &[Complex],
    zeros: &[Complex],
    a0: f32,
    power: i32,
) -> Vec<Complex> {
    let mut z = vec![Complex::default()];
    for i in 1..nt / 2 {
        let w = dw * i as f32;
        let mut zzz = c(1.0);
        for &pole in poles {
            zzz = zzz / imag_minus(w, pole);
        }
        for &zero in zeros {
            zzz = zzz * imag_minus(w, zero);
        }
        z.push((c(a0) * zzz) * Complex::new(0.0, w).powi(power));
    }
    z
}

/// The WWSSN long-period response for a delta input, normalised to a peak
/// gain of one. `ts`, `tg` are the seismometer and galvanometer periods,
/// `hs`, `hg` their damping, and `s` the coupling.
pub fn instw(nt: usize, dw: f32, ts: f32, tg: f32, hs: f32, hg: f32, s: f32) -> Vec<Complex> {
    let ws = (2.0 * PI) / ts;
    let wg = (2.0 * PI) / tg;
    let a0 = (ws * wg) * (ws * wg);
    let a1 = ((2.0 * ws) * wg) * (hs * wg + hg * ws);
    let a2 = (ws * ws + wg * wg) + ((((4.0 * (1.0 - s)) * ws) * wg) * hs) * hg;
    let a3 = 2.0 * (hs * ws + hg * wg);
    let den = |zw: Complex| c(a0) + zw * (c(a1) + zw * (c(a2) + zw * (c(a3) + zw)));
    let omegas = || (1..nt / 2).map(|i| Complex::new(0.0, dw * i as f32));
    let mut gm = 0.0f32;
    for zw in omegas() {
        gm = gm.max((((zw * zw) * zw) / den(zw)).abs());
    }
    let g = 1.0 / gm;
    let mut z = vec![Complex::default()];
    // The reference raises to `3 - id` with `id` a dummy argument, so GCC
    // calls `pow_c4_i4` instead of multiplying out a constant power.
    z.extend(omegas().map(|zw| (zw.powi(3) * c(g)) / den(zw)));
    z
}
