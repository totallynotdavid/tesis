//! Haskell propagator matrices: the response of a layered structure near the
//! source (`refl`, `reflsh`) and near the receiver (`cnvr`, `cnvrsh`) to a
//! plane wave of ray parameter `p`, per frequency of a half spectrum.
//!
//! The last layer of each structure is the half-space. A layer with `vs = 0`
//! is water.

use super::c;
use crate::fft::Complex;
use crate::formats::structure::Layer;

type Matrix<const N: usize> = [[Complex; N]; N];

/// Near-source coefficients of the up- and down-going P and S waves in the
/// source layer, relative to the wave leaving the half-space. P-SV fills
/// all four; SH fills `su` and `sd`.
pub struct Source {
    pub pu: Vec<Complex>,
    pub pd: Vec<Complex>,
    pub su: Vec<Complex>,
    pub sd: Vec<Complex>,
    /// Vertical travel time from the top of the source layer to the
    /// half-space: P for P-SV, S for SH.
    pub time: f32,
}

/// Free-surface displacement per unit incident displacement from the
/// half-space.
pub struct Receiver {
    pub horizontal: Vec<Complex>,
    /// Upward.
    pub vertical: Vec<Complex>,
    /// Vertical travel time from the half-space to the surface.
    pub time: f32,
}

/// P-SV near-source coefficients for source layer `l`, for an incident P
/// (`p_wave`) or SV wave.
pub fn refl(layers: &[Layer], l: usize, p_wave: bool, p: f32, dw: f32, half: usize) -> Source {
    let nl = layers.len();
    let p2 = p * p;
    let time = p_time(&layers[l..nl - 1], p2);
    let zel = inverse_e(&layers[l], p, p2);
    let ze = inverse_e_half_space(&layers[nl - 1], p, p2);
    let mut out = Source {
        pu: Vec::with_capacity(half),
        pd: Vec::with_capacity(half),
        su: Vec::with_capacity(half),
        sd: Vec::with_capacity(half),
        time,
    };
    for i in 0..half {
        let w = i as f32 * dw;
        let mut zaa = identity::<4>();
        let mut zjl = if l == 0 {
            prod(&zel, &zaa)
        } else {
            Matrix::default()
        };
        for (m, layer) in layers[..nl - 1].iter().enumerate() {
            zaa = prod(&layer_matrix(layer, p, p2, w), &zaa);
            if m + 1 == l {
                zjl = prod(&zel, &zaa);
            }
        }
        let zj = prod(&ze, &zaa);
        let [zj1, zj2, zj3, zj4] = j_minors(&zj);
        let zdet = zj4 * zj1 - zj3 * zj2;
        let l = &zjl;
        let two = c(2.0);
        if p_wave {
            out.pu
                .push((zj1 * (l[0][0] + l[1][0]) + zj2 * (l[0][1] + l[1][1])) / zdet);
            out.pd
                .push((zj1 * (l[0][0] - l[1][0]) + zj2 * (l[0][1] - l[1][1])) / zdet);
            out.su
                .push((two * (zj1 * (l[2][0] + l[3][0]) + zj2 * (l[2][1] + l[3][1]))) / zdet);
            out.sd
                .push((two * (zj1 * (l[3][0] - l[2][0]) + zj2 * (l[3][1] - l[2][1]))) / zdet);
        } else {
            out.pu
                .push((zj3 * (l[0][0] + l[1][0]) + zj4 * (l[0][1] + l[1][1])) / zdet / two);
            out.pd
                .push((zj3 * (l[0][0] - l[1][0]) + zj4 * (l[0][1] - l[1][1])) / zdet / two);
            out.su
                .push((zj3 * (l[2][0] + l[3][0]) + zj4 * (l[2][1] + l[3][1])) / zdet);
            out.sd
                .push((zj3 * (l[3][0] - l[2][0]) + zj4 * (l[3][1] - l[2][1])) / zdet);
        }
    }
    out
}

/// P-SV free-surface conversion for an incident P (`p_wave`) or SV wave.
pub fn cnvr(layers: &[Layer], p_wave: bool, p: f32, dw: f32, half: usize) -> Receiver {
    let nl = layers.len();
    let p2 = p * p;
    let bottom = &layers[nl - 1];
    let ze = inverse_e_half_space(bottom, p, p2);
    let mut out = Receiver {
        horizontal: Vec::with_capacity(half),
        vertical: Vec::with_capacity(half),
        time: p_time(&layers[..nl - 1], p2),
    };
    for i in 0..half {
        let w = i as f32 * dw;
        let mut zaa = identity::<4>();
        for layer in &layers[..nl - 1] {
            zaa = prod(&layer_matrix(layer, p, p2, w), &zaa);
        }
        let zj = prod(&ze, &zaa);
        let [zj1, zj2, zj3, zj4] = j_minors(&zj);
        let zdet = zj4 * zj1 - zj3 * zj2;
        if p_wave {
            let d = (zdet * c(bottom.vp)) * c(p);
            out.horizontal.push(-((zj1 * c(2.0)) / d));
            out.vertical.push((zj2 * c(2.0)) / d);
        } else {
            let d = (zdet * c(bottom.vs)) * c(p);
            out.horizontal.push(zj3 / d);
            out.vertical.push(zj4 / d);
        }
    }
    out
}

/// SH near-source coefficients for source layer `l`.
pub fn reflsh(layers: &[Layer], l: usize, p: f32, dw: f32, half: usize) -> Source {
    let nl = layers.len();
    let p2 = p * p;
    let rgb = sh_rigidity(&layers[nl - 1], p, p2);
    let rgbl = sh_rigidity(&layers[l], p, p2);
    let mut out = Source {
        pu: Vec::new(),
        pd: Vec::new(),
        su: Vec::with_capacity(half),
        sd: Vec::with_capacity(half),
        time: s_time(&layers[l..nl - 1], p2),
    };
    for i in 0..half {
        let w = i as f32 * dw;
        let mut zaa = identity::<2>();
        let mut zal = if l == 0 {
            identity::<2>()
        } else {
            Matrix::default()
        };
        for (m, layer) in layers[..nl - 1].iter().enumerate() {
            let za = if layer.vs != 0.0 {
                sh_matrix(layer, p, p2, w)
            } else {
                [[c(1.0), c(0.0)], [c(0.0), c(0.0)]]
            };
            zaa = prod(&za, &zaa);
            if m + 1 == l {
                zal = zaa;
            }
        }
        let zdet = zaa[0][0] + zaa[1][0] / c(rgb);
        out.su.push((zal[0][0] - zal[1][0] / c(rgbl)) / zdet);
        out.sd.push((zal[0][0] + zal[1][0] / c(rgbl)) / zdet);
    }
    out
}

/// SH free-surface conversion; the response is horizontal.
pub fn cnvrsh(layers: &[Layer], p: f32, dw: f32, half: usize) -> Receiver {
    let nl = layers.len();
    let p2 = p * p;
    let rgb = sh_rigidity(&layers[nl - 1], p, p2);
    let mut out = Receiver {
        horizontal: Vec::with_capacity(half),
        vertical: Vec::new(),
        time: s_time(&layers[..nl - 1], p2),
    };
    for i in 0..half {
        let w = i as f32 * dw;
        let mut zaa = identity::<2>();
        // The reference does not treat water here.
        for layer in &layers[..nl - 1] {
            zaa = prod(&sh_matrix(layer, p, p2, w), &zaa);
        }
        let zdet = zaa[0][0] + zaa[1][0] / c(rgb);
        out.horizontal.push(c(2.0) / zdet);
    }
    out
}

fn p_time(layers: &[Layer], p2: f32) -> f32 {
    let mut t = 0.0;
    for layer in layers {
        t += layer.thickness * (1.0 / (layer.vp * layer.vp) - p2).sqrt();
    }
    t
}

fn s_time(layers: &[Layer], p2: f32) -> f32 {
    let mut t = 0.0;
    for layer in layers.iter().filter(|layer| layer.vs != 0.0) {
        t += layer.thickness * (1.0 / (layer.vs * layer.vs) - p2).sqrt();
    }
    t
}

fn identity<const N: usize>() -> Matrix<N> {
    let mut m = [[Complex::default(); N]; N];
    for (j, row) in m.iter_mut().enumerate() {
        row[j] = c(1.0);
    }
    m
}

/// `a b`, each sum accumulated from zero as `sub.prod.f` does.
fn prod<const N: usize>(a: &Matrix<N>, b: &Matrix<N>) -> Matrix<N> {
    let mut out = [[Complex::default(); N]; N];
    for (j1, row) in out.iter_mut().enumerate() {
        for (j2, cell) in row.iter_mut().enumerate() {
            let mut s = Complex::default();
            for j in 0..N {
                s = s + a[j1][j] * b[j][j2];
            }
            *cell = s;
        }
    }
    out
}

fn j_minors(zj: &Matrix<4>) -> [Complex; 4] {
    [
        zj[3][1] - zj[2][1],
        zj[2][0] - zj[3][0],
        zj[1][1] - zj[0][1],
        zj[0][0] - zj[1][0],
    ]
}

/// The P-SV layer matrix of `layer` at angular frequency `w`.
fn layer_matrix(layer: &Layer, p: f32, p2: f32, w: f32) -> Matrix<4> {
    let vp2 = 1.0 / (layer.vp * layer.vp);
    let y1 = (vp2 - p2).sqrt();
    let pm = (y1 * w) * layer.thickness;
    let (spm, cpm) = pm.sin_cos();
    let zspm = Complex::new(0.0, spm);
    let ra = y1 / p;
    let rmc2 = layer.density / p2;
    if layer.vs == 0.0 {
        let zero = Complex::default();
        let a12 = -(zspm / c(ra));
        let a13 = c(-((1.0 / rmc2) * cpm));
        let a14 = (c(1.0 / rmc2) * zspm) / c(ra);
        return [
            [c(1.0), a12, a13, a14],
            [zero, c(cpm), c(ra / rmc2) * zspm, a13],
            [zero, c(rmc2 / ra) * zspm, c(cpm), a12],
            [zero; 4],
        ];
    }
    let vs2 = 1.0 / (layer.vs * layer.vs);
    let y2 = (vs2 - p2).sqrt();
    let gm = 2.0 * ((p * layer.vs) * (p * layer.vs));
    let qm = (y2 * w) * layer.thickness;
    let (sqm, cqm) = qm.sin_cos();
    let zsqm = Complex::new(0.0, sqm);
    let rb = y2 / p;
    let g1 = gm - 1.0;
    let a11 = c(gm * cpm - g1 * cqm);
    let a12 = c(g1 / ra) * zspm + c(gm * rb) * zsqm;
    let a13 = c(-((1.0 / rmc2) * (cpm - cqm)));
    let a14 = c(1.0 / rmc2) * (zspm / c(ra) + c(rb) * zsqm);
    let a21 = -(c(gm * ra) * zspm + c(g1 / rb) * zsqm);
    let a22 = c(-(g1 * cpm) + gm * cqm);
    let a23 = c(1.0 / rmc2) * (c(ra) * zspm + zsqm / c(rb));
    let a31 = c(((rmc2 * gm) * g1) * (cpm - cqm));
    let a32 = c(rmc2) * (c((g1 * g1) / ra) * zspm + c((gm * gm) * rb) * zsqm);
    let a41 = c(rmc2) * (c((gm * gm) * ra) * zspm + c((g1 * g1) / rb) * zsqm);
    [
        [a11, a12, a13, a14],
        [a21, a22, a23, a13],
        [a31, a32, a22, a12],
        [a41, a31, a21, a11],
    ]
}

/// E⁻¹ of the source layer, which takes the layer's stress-displacement
/// vector to its wave potentials.
fn inverse_e(layer: &Layer, p: f32, p2: f32) -> Matrix<4> {
    let (vp, vs, den) = (layer.vp, layer.vs, layer.density);
    let gm = 2.0 * ((p * vs) * (p * vs));
    let ra = (1.0 / (vp * vp) - p2).sqrt() / p;
    let rb = (1.0 / (vs * vs) - p2).sqrt() / p;
    let mut e = inverse_e_common(layer, gm, rb, p2);
    e[1][1] = c((gm - 1.0) / (((vp * vp) * ra) * p2));
    e[1][3] = c(1.0 / ((den * (vp * vp)) * ra));
    e
}

/// E⁻¹ of the half-space, where P may be evanescent: past the P critical
/// slowness its vertical slowness is imaginary.
fn inverse_e_half_space(layer: &Layer, p: f32, p2: f32) -> Matrix<4> {
    let (vp, vs, den) = (layer.vp, layer.vs, layer.density);
    let gm = 2.0 * ((p * vs) * (p * vs));
    let vp2 = 1.0 / (vp * vp);
    let zy1 = if vp2 >= p2 {
        c((vp2 - p2).sqrt())
    } else {
        Complex::new(0.0, -(p2 - vp2).sqrt())
    };
    let zra = zy1 / c(p);
    let rb = (1.0 / (vs * vs) - p2).sqrt() / p;
    let mut e = inverse_e_common(layer, gm, rb, p2);
    e[1][1] = c(gm - 1.0) / ((c(vp * vp) * zra) * c(p2));
    e[1][3] = c(1.0) / (c(den * (vp * vp)) * zra);
    e
}

fn inverse_e_common(layer: &Layer, gm: f32, rb: f32, p2: f32) -> Matrix<4> {
    let (vp, vs, den) = (layer.vp, layer.vs, layer.density);
    let mut e = Matrix::<4>::default();
    e[0][0] = c(-(2.0 * ((vs / vp) * (vs / vp))));
    e[0][2] = c(1.0 / (den * (vp * vp)));
    e[2][0] = c((gm - 1.0) / (gm * rb));
    e[2][2] = c(-(p2 / ((den * gm) * rb)));
    e[3][1] = c(1.0);
    e[3][3] = c(p2 / (den * gm));
    e
}

/// The SH impedance term `den vs² ηβ / p`.
fn sh_rigidity(layer: &Layer, p: f32, p2: f32) -> f32 {
    let y2 = (1.0 / (layer.vs * layer.vs) - p2).sqrt();
    ((layer.density * (layer.vs * layer.vs)) * y2) / p
}

fn sh_matrix(layer: &Layer, p: f32, p2: f32, w: f32) -> Matrix<2> {
    let y2 = (1.0 / (layer.vs * layer.vs) - p2).sqrt();
    let rgb = ((layer.density * (layer.vs * layer.vs)) * y2) / p;
    let qm = (y2 * w) * layer.thickness;
    let (sqm, cqm) = qm.sin_cos();
    let zsqm = Complex::new(0.0, sqm);
    [[c(cqm), c(1.0 / rgb) * zsqm], [c(rgb) * zsqm, c(cqm)]]
}
