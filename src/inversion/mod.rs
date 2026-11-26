//! `inversion`: the multiple-shock grid search. Each subevent is the grid
//! point, onset and moment tensor whose synthetics best fit the residual
//! records. With more than one subevent, two more rounds refit each one
//! against the records with the others removed.
//!
//! The port matches the gfortran build bit for bit, so it keeps the
//! reference's f32 constants, its order of operations, and comparisons
//! written so that NaN takes the same branch.
#![allow(
    clippy::approx_constant,
    clippy::needless_range_loop,
    clippy::neg_cmp_op_on_partial_ord,
    clippy::assign_op_pattern,
    reason = "bit parity with the reference build"
)]

pub mod linalg;
pub mod mechanism;

use std::fmt;

use crate::formats::greens::Greens;
use crate::formats::inversion::{Basis, Parameters, Shape, TimeFunction};
use crate::formats::solution::{Fit, Mechanism, Solution, Subevent};
use crate::formats::waveforms::{Instrument, Station, Waveforms};
use crate::green;
use linalg::{Matrix, N};
use mechanism::{DoubleCouple, Weights, Zmom};

const RAD: f32 = 0.01745329;

#[derive(Debug, PartialEq)]
pub struct Error {
    pub message: String,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for Error {}

fn error(message: String) -> Error {
    Error { message }
}

/// A station the inversion fits.
struct Record<'a> {
    station: &'a Station,
    instrument: i32,
    weight: f32,
    onset: i32,
    velocity: f32,
    /// The record at the common gain, times the weight.
    observed: Vec<f32>,
    /// The convolved base elements, `[depth][element]`, at the same gain.
    elements: Vec<[Vec<f32>; N]>,
    /// The delay of each grid point, `[point][depth]`, in samples.
    delays: Vec<Vec<i32>>,
}

/// The sizes and constants every stage shares.
struct Setup {
    nt: usize,
    nm: usize,
    dt: f32,
    points: Vec<(f32, f32)>,
}

pub fn invert(p: &Parameters, waveforms: &Waveforms, greens: &Greens) -> Result<Solution, Error> {
    let dt = greens.dt;
    let nt = (p.length / dt + 0.5) as i32;
    if nt < 1 {
        return Err(error(format!("length {} s is under one sample", p.length)));
    }
    let setup = Setup {
        nt: nt as usize,
        nm: p.basis.elements(),
        dt,
        points: p.grid.points(),
    };
    let source = stime(setup.nt, dt, p.time_function);
    let records = records(p, waveforms, greens, &setup, &source)?;
    let mut rww = vec![[[0f32; N]; N]; greens.nk];
    let mut xpower = 0f32;
    for r in &records {
        for (k, elements) in r.elements.iter().enumerate() {
            for n1 in 0..setup.nm {
                for n2 in 0..setup.nm {
                    for i in 0..setup.nt {
                        rww[k][n1][n2] += elements[n1][i] * elements[n2][i];
                    }
                }
            }
        }
        for x in &r.observed {
            xpower += x * x;
        }
    }
    let arww = rww
        .iter()
        .enumerate()
        .map(|(k, m)| {
            inverse(m, setup.nm)
                .ok_or_else(|| error(format!("depth {}: singular normal equations", k + 1)))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let search = search(p, greens, &setup, &records, &rww, &arww, xpower)?;
    solution(p, greens, &setup, source, &records, &arww, search, xpower)
}

/// The stations `green` kept, in order, less those with weight 0.
fn records<'a>(
    p: &Parameters,
    waveforms: &'a Waveforms,
    greens: &Greens,
    setup: &Setup,
    source: &[f32],
) -> Result<Vec<Record<'a>>, Error> {
    let mut traces = greens.stations.iter();
    let mut weights = p.weights.iter().copied().chain(std::iter::repeat(0.0));
    let mut records = Vec::new();
    for station in &waveforms.stations {
        if green::is_skipped(station) {
            continue;
        }
        let Some(trace) = traces.next() else { break };
        let weight = weights.next().unwrap_or(0.0);
        if weight == 0.0 {
            continue;
        }
        let name = || green::short_name(station);
        let (instrument, gain) = match station.instrument {
            Instrument::Wwssn { magnification, .. } => (0, magnification),
            Instrument::Gdsn { magnification, .. } => (1, magnification * 1.0e4),
        };
        let wt = (1.0e4 / gain) * weight;
        let onset = station.onset + p.shift;
        let observed = (0..setup.nt)
            .map(|i| {
                let i1 = ((onset as f32 + (i as f32 * setup.dt) / station.dt) + 1.5) as i32;
                if i1 > 0 && i1 as usize <= station.samples.len() {
                    (station.samples[i1 as usize - 1] * station.peak) * wt
                } else {
                    0.0
                }
            })
            .collect();
        let elements = trace
            .depths
            .iter()
            .map(|depth| {
                std::array::from_fn(|n| {
                    if n >= setup.nm {
                        return Vec::new();
                    }
                    conv(&depth[n], source, setup.nt)
                        .into_iter()
                        .map(|y| (y * wt) * setup.dt)
                        .collect()
                })
            })
            .collect();
        let velocity = match station.phase {
            1 | 4 => greens.vp,
            2 | 3 => greens.vs,
            phase => {
                return Err(error(format!(
                    "station {}: phase {phase} is not 1 to 4",
                    name().trim_end()
                )));
            }
        };
        let delays = setup
            .points
            .iter()
            .map(|&(fr, ffi)| {
                (1..=greens.nk as i32)
                    .map(|k| delay(station, velocity, greens, fr, ffi, k) / setup.dt)
                    .map(|d| d.round() as i32)
                    .collect()
            })
            .collect();
        records.push(Record {
            station,
            instrument,
            weight,
            onset,
            velocity,
            observed,
            elements,
            delays,
        });
    }
    Ok(records)
}

/// The travel-time difference, in seconds, between grid point `(fr, ffi)`
/// at depth `k` and the hypocentre.
fn delay(station: &Station, v: f32, greens: &Greens, fr: f32, ffi: f32, k: i32) -> f32 {
    let p = station.slowness;
    let az = (station.azimuth - ffi) * RAD;
    let dip = greens.dip * RAD;
    -((fr * p) * az.cos())
        - (greens.dk * (k - greens.k0) as f32)
            * (((1.0 / (v * v) - p * p).sqrt() * dip.sin()) + (p * dip.cos()) * az.sin())
}

/// The inverse of the leading `nm` block, zero elsewhere.
fn inverse(m: &Matrix, nm: usize) -> Option<Matrix> {
    let mut lu = *m;
    let indx = linalg::ludcmp(&mut lu, nm)?;
    let mut inverse = [[0f32; N]; N];
    for n in 0..nm {
        let mut b = [0f32; N];
        b[n] = 1.0;
        linalg::lubksb(&lu, nm, &indx, &mut b);
        for n1 in 0..nm {
            inverse[n1][n] = b[n1];
        }
    }
    Some(inverse)
}

/// A subevent's grid point: onset sample from 1, point and depth from 0.
#[derive(Clone, Copy)]
struct Best {
    i: i32,
    l: usize,
    k: usize,
}

struct Search {
    rounds: Vec<Vec<Subevent>>,
    weights: Vec<Weights>,
    convergence: Vec<f32>,
    /// The residuals, as `observed`.
    residuals: Vec<Vec<f32>>,
    xp: f32,
}

fn search(
    p: &Parameters,
    greens: &Greens,
    setup: &Setup,
    records: &[Record],
    rww: &[Matrix],
    arww: &[Matrix],
    xpower: f32,
) -> Result<Search, Error> {
    let (nt, nm, dt) = (setup.nt, setup.nm, setup.dt);
    let ne = p.spans.len();
    let rounds = if ne > 1 { 3 } else { 1 };
    let mut x: Vec<Vec<f32>> = records.iter().map(|r| r.observed.clone()).collect();
    let mut syn = vec![vec![vec![0f32; nt]; records.len()]; ne];
    let mut rwx = vec![vec![[const { Vec::new() }; N]; greens.nk]; records.len()];
    let mut aop = vec![[0f32; N]; ne];
    let mut best: Vec<Option<Best>> = vec![None; ne];
    let mut eror = vec![1f32];
    // The reference keeps these between grid points and subevents.
    let mut zmom = Zmom::default();
    let mut il = 0;
    let mut rh1 = [0f32; N];
    let mut xp = 0f32;
    let mut out = Vec::new();
    for round in 1..=rounds {
        let mut subevents = Vec::new();
        for ns in 0..ne {
            let mut opt = 0f32;
            let (mut rww0, mut rwx0) = (0f32, 0f32);
            for (j, r) in records.iter().enumerate() {
                if round != 1 {
                    for i in 0..nt {
                        let s = syn[ns][j][i];
                        x[j][i] += s;
                        rww0 += s * s;
                        rwx0 += x[j][i] * s;
                    }
                }
                for (k, elements) in r.elements.iter().enumerate() {
                    for n in 0..nm {
                        rwx[j][k][n] = corr(&elements[n], &x[j]);
                    }
                }
            }
            if round >= 3 {
                opt = (rwx0 * rwx0) / rww0;
                let am = rwx0 / rww0;
                for v in &mut aop[ns][..nm] {
                    *v = am * *v;
                }
            }
            for k in 0..greens.nk {
                let kk = k as i32 + 1 - greens.k0;
                for (l, &(fr, _)) in setup.points.iter().enumerate() {
                    let is1 = ((((fr * fr + (kk * kk) as f32 * (greens.dk * greens.dk)).sqrt()
                        / p.max_velocity)
                        / dt)
                        + 1.5) as i32;
                    let is2 = (p.spans[ns] / dt + 1.5) as i32;
                    for i in is1..=is2.min(nt as i32) {
                        let mut rh0 = [0f32; N];
                        for (n, sum) in rh0.iter_mut().enumerate().take(nm) {
                            for (j, r) in records.iter().enumerate() {
                                let at = i + r.delays[l][k];
                                if at >= 1 && at as usize <= nt {
                                    *sum += rwx[j][k][n][at as usize - 1];
                                }
                            }
                        }
                        let mut rh = [0f32; N];
                        for n in 0..nm {
                            for n1 in 0..nm {
                                rh[n] += arww[k][n][n1] * rh0[n1];
                            }
                        }
                        let mut obj = 0f32;
                        for n1 in 0..nm {
                            obj += rh0[n1] * rh[n1];
                        }
                        rh1[..nm].copy_from_slice(&rh[..nm]);
                        if p.basis.is_double_couple() {
                            if let Some(flag) = zmom.solve(&arww[k], &rww[k], &rh, &mut rh1) {
                                il = flag;
                            }
                            if il != 0 {
                                continue;
                            }
                            let mut dobj = 0f32;
                            for n1 in 0..nm {
                                for n2 in 0..nm {
                                    dobj +=
                                        (rww[k][n1][n2] * (rh1[n1] - rh[n1])) * (rh1[n2] - rh[n2]);
                                }
                            }
                            obj -= dobj;
                            if obj < 0.0 {
                                obj = 0.0;
                            }
                        }
                        if !(obj <= opt) {
                            best[ns] = Some(Best { i, l, k });
                            opt = obj;
                            aop[ns][..nm].copy_from_slice(&rh1[..nm]);
                        }
                    }
                }
            }
            let b = best[ns]
                .ok_or_else(|| error(format!("subevent {}: no grid point fits", ns + 1)))?;
            if round == 1 {
                eror.push(eror[ns] - opt / xpower);
            }
            let dc = dcp(&aop[ns])?;
            xp = 0.0;
            for (j, r) in records.iter().enumerate() {
                let isf = b.i + r.delays[b.l][b.k] - 1;
                for i in 1..=nt as i32 {
                    let i1 = i - isf;
                    let mut s = 0f32;
                    if i1 >= 1 && i1 as usize <= nt {
                        for n1 in 0..nm {
                            s += aop[ns][n1] * r.elements[b.k][n1][i1 as usize - 1];
                        }
                        x[j][i as usize - 1] -= s;
                    }
                    syn[ns][j][i as usize - 1] = s;
                    let v = x[j][i as usize - 1];
                    xp += v * v;
                }
            }
            let (distance, azimuth) = setup.points[b.l];
            subevents.push(Subevent {
                time: dt * (b.i - 1) as f32,
                distance,
                azimuth,
                offset: greens.dk * (b.k as i32 + 1 - greens.k0) as f32,
                mechanism: mechanism(&dc, &aop[ns], p.basis),
                variance: xp / xpower,
                tensor: mechanism::columns(&aop[ns]),
                axis: dc.axis,
            });
        }
        out.push(subevents);
    }
    Ok(Search {
        rounds: out,
        weights: aop,
        convergence: eror[1..].to_vec(),
        residuals: x,
        xp,
    })
}

#[allow(clippy::too_many_arguments)]
fn solution(
    p: &Parameters,
    greens: &Greens,
    setup: &Setup,
    source: Vec<f32>,
    records: &[Record],
    arww: &[Matrix],
    search: Search,
    xpower: f32,
) -> Result<Solution, Error> {
    let nm = setup.nm;
    let mut total = [0f32; N];
    for weights in &search.weights {
        for n in 0..nm {
            total[n] += weights[n];
        }
    }
    let nfree = setup.nt as i32 * records.len() as i32 - (nm * p.spans.len()) as i32;
    let sx = search.xp / nfree as f32;
    let at_h0 = usize::try_from(greens.k0 - 1)
        .ok()
        .and_then(|k| arww.get(k))
        .copied()
        .unwrap_or_default();
    let mut total1 = [0f32; N];
    for n in 0..nm {
        total1[n] = total[n] + (at_h0[n][n] * sx).sqrt().abs().copysign(total[n]);
    }
    let dc1 = dcp(&total1)?;
    let dc = dcp(&total)?;
    let m1 = mechanism(&dc1, &total1, p.basis);
    let m = mechanism(&dc, &total, p.basis);
    let stations = records
        .iter()
        .zip(&search.residuals)
        .map(|(r, residual)| Fit {
            name: green::short_name(r.station),
            instrument: r.instrument,
            phase: r.station.phase,
            component: r.station.component,
            azimuth: r.station.azimuth,
            back_azimuth: r.station.back_azimuth,
            distance: r.station.distance,
            slowness: r.station.slowness,
            spreading: r.station.spreading,
            weight: r.weight,
            onset: r.onset,
            velocity: r.velocity,
            observed: r.observed.iter().map(|x| x / r.weight).collect(),
            residual: residual.iter().map(|x| x / r.weight).collect(),
        })
        .collect();
    Ok(Solution {
        title: p.title.clone(),
        structure: greens.title.clone(),
        dt: setup.dt,
        h0: greens.h0,
        dip: greens.dip,
        basis: p.basis,
        shift: p.shift,
        time_function: p.time_function,
        max_velocity: p.max_velocity,
        spans: p.spans.clone(),
        source,
        dk: greens.dk,
        nk: greens.nk,
        k0: greens.k0,
        points: setup.points.clone(),
        stations,
        rounds: search.rounds,
        convergence: search.convergence,
        variance: search.xp / xpower,
        total: m,
        deviation: Mechanism {
            moment: m1.moment - m.moment,
            second: m1.second - m.second,
            strike: m1.strike - m.strike,
            dip: m1.dip - m.dip,
            slip: m1.slip - m.slip,
        },
        tensor: mechanism::columns(&total),
        axis: dc.axis,
        magnitude: (m.moment.log10() + 8.9) / 1.5,
    })
}

fn dcp(v: &Weights) -> Result<DoubleCouple, Error> {
    mechanism::dcp(v).ok_or_else(|| error("eigenvalues of a moment tensor did not converge".into()))
}

/// The mechanism a report prints. Under [`Basis::MomentTensor`] the second
/// value is the isotropic weight rather than the non-double-couple part.
fn mechanism(dc: &DoubleCouple, v: &Weights, basis: Basis) -> Mechanism {
    Mechanism {
        moment: dc.moment,
        second: if basis == Basis::MomentTensor {
            v[5]
        } else {
            dc.non_dc
        },
        strike: dc.strike,
        dip: dc.dip,
        slip: dc.slip,
    }
}

/// `stime`: the source time function, `n` samples of `dt`. A flank longer
/// than `n` samples is cut at `n`.
fn stime(n: usize, dt: f32, t: TimeFunction) -> Vec<f32> {
    const PI: f32 = 3.1415926;
    let (t1, t2) = (t.t1, t.t2);
    let mut x = vec![0f32; n];
    let it1 = (t1 / dt + 1.5) as i32;
    let it2 = (t2 / dt + 1.5) as i32;
    let it3 = ((t1 + t2) / dt + 1.5) as i32;
    // Samples `from..=to`, counted from 1, as a function of `i - 1`.
    let mut fill = |from: i32, to: i32, value: &dyn Fn(f32) -> f32| {
        for i in from.max(1)..=to.min(n as i32) {
            x[i as usize - 1] = value((i - 1) as f32);
        }
    };
    let n = n as i32;
    match t.shape {
        Shape::Impulse => fill(1, 1, &|_| 1.0 / dt),
        Shape::Trapezoid if t2 < t1 => {
            fill(2, it1, &|m| (m * dt) / t1);
            fill(it1 + 1, n, &|_| 1.0);
        }
        Shape::Trapezoid => {
            fill(2, it1, &|m| (m * dt) / (t1 * t2));
            fill(it1 + 1, it2, &|_| 1.0 / t2);
            fill(it2 + 1, it3, &|m| ((t1 + t2) - m * dt) / (t1 * t2));
        }
        Shape::Cosine if t2 < t1 => {
            fill(2, it1, &|m| (1.0 - (((PI * m) * dt) / t1).cos()) * 0.5);
            fill(it1 + 1, n, &|_| 1.0);
        }
        Shape::Cosine => {
            fill(2, it1, &|m| {
                (1.0 - (((PI * m) * dt) / t1).cos()) / (2.0 * t2)
            });
            fill(it1 + 1, it2, &|_| 1.0 / t2);
            fill(it2 + 1, it3, &|m| {
                (1.0 - ((PI * ((t1 + t2) - m * dt)) / t1).cos()) / (2.0 * t2)
            });
        }
        Shape::Triangle => {
            fill(2, it1, &|m| ((2.0 * m) * dt) / (t1 * t2));
            fill(it1 + 1, it2, &|m| ((2.0 * (t2 - m * dt)) / (t2 - t1)) / t2);
        }
    }
    x
}

/// `conv`: the first `n` samples of `x * y`.
fn conv(x: &[f32], y: &[f32], n: usize) -> Vec<f32> {
    (0..n)
        .map(|i| {
            let mut z = 0f32;
            for j in 0..x.len().min(i + 1) {
                if let Some(&v) = y.get(i - j) {
                    z += x[j] * v;
                }
            }
            z
        })
        .collect()
}

/// `corr`: the correlation of `x` with `y` at lags `0..y.len()`, both the
/// same length.
fn corr(x: &[f32], y: &[f32]) -> Vec<f32> {
    (0..y.len())
        .map(|i| {
            let mut z = 0f32;
            for (j, &a) in x.iter().enumerate() {
                if let Some(&b) = y.get(i + j) {
                    z += a * b;
                }
            }
            z
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_trapezoid_has_unit_area() {
        let t = TimeFunction {
            shape: Shape::Trapezoid,
            t1: 6.0,
            t2: 10.0,
        };
        let area: f32 = stime(160, 0.5, t).iter().sum::<f32>() * 0.5;
        assert!((area - 1.0).abs() < 1e-5, "{area}");
    }

    #[test]
    fn conv_and_corr_shift_by_one_sample() {
        assert_eq!(conv(&[0.0, 1.0], &[1.0, 2.0, 3.0], 3), [0.0, 1.0, 2.0]);
        assert_eq!(corr(&[0.0, 1.0, 0.0], &[1.0, 2.0, 3.0]), [2.0, 3.0, 0.0]);
    }
}
