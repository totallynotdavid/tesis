//! Green's functions for teleseismic body waves: for each station and each
//! point of a line of sources, the synthetic seismogram of a set of
//! mechanisms. `green.f` uses the six base moment-tensor elements, `green3.f`
//! one fixed mechanism.
//!
//! A synthetic combines the radiation pattern, the near-source structure
//! with its depth phases, the near-receiver structure, Q and the
//! instrument, in the frequency domain, and transforms it back to time.

mod propagator;
mod radiation;
mod response;

use std::fmt;

use crate::fft::{Complex, PI, cfft};
use crate::formats::fixed;
use crate::formats::greens::{self, ELEMENTS};
use crate::formats::structure::{Layer, Structure};
use crate::formats::waveforms::{Instrument, Station, Waveforms};
use propagator::{Receiver, Source};
use radiation::{Radiation, radp};

/// The reference's degrees-to-radians factor.
pub(crate) const RAD: f32 = 0.0174533;

/// A double couple, in degrees. A dip above 360 selects the isotropic
/// source.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mechanism {
    pub strike: f32,
    pub dip: f32,
    pub rake: f32,
}

/// The six base elements in `fort.3` order: Mxy, Mzz-Myy, Myz, Mzx, Mzz-Mxx
/// and the isotropic Mxx+Myy+Mzz.
pub const BASE_ELEMENTS: [Mechanism; ELEMENTS] = [
    Mechanism::new(0.0, 90.0, 0.0),
    Mechanism::new(135.0, 90.0, 0.0),
    Mechanism::new(180.0, 90.0, 90.0),
    Mechanism::new(90.0, 90.0, 90.0),
    Mechanism::new(90.0, 45.0, 90.0),
    Mechanism::new(0.0, 500.0, 0.0),
];

impl Mechanism {
    pub const fn new(strike: f32, dip: f32, rake: f32) -> Self {
        Self { strike, dip, rake }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sampling {
    /// Samples per trace, a power of two.
    pub nt: usize,
    pub dt: f32,
}

/// Source points on a line through the reference point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Line {
    /// The reference depth in km, at point `k0`, counted from 1.
    pub h0: f32,
    /// Spacing of the points along the line, in km.
    pub dk: f32,
    pub nk: usize,
    pub k0: i32,
    /// Dip of the line in degrees.
    pub dip: f32,
}

impl Line {
    /// The depth of point `k`, counted from 1.
    fn depth(&self, k: usize) -> f32 {
        self.h0 + (((k as i32 - self.k0) as f32 * self.dk) * (self.dip * RAD).sin())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Phase {
    P,
    Sv,
    Sh,
    /// P reflected once at the bounce point. It arrives 10 s late and is
    /// Hilbert transformed.
    Pp,
}

impl Phase {
    fn from_code(ib: i32) -> Option<Self> {
        match ib {
            1 => Some(Self::P),
            2 => Some(Self::Sv),
            3 => Some(Self::Sh),
            4 => Some(Self::Pp),
            _ => None,
        }
    }

    fn is_p(self) -> bool {
        matches!(self, Self::P | Self::Pp)
    }
}

#[derive(Debug, PartialEq)]
pub struct Error {
    pub station: String,
    pub message: String,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "station {}: {}", self.station, self.message)
    }
}

impl std::error::Error for Error {}

/// `green` leaves out a station whose onset is past its last sample.
pub fn is_skipped(station: &Station) -> bool {
    station.onset > station.samples.len() as i32
}

/// The name `fort.3` keeps for a station: the first 12 characters.
pub fn short_name(station: &Station) -> String {
    station.name.chars().take(12).collect()
}

/// `green`: the base elements, for `fort.3`.
pub fn base(
    structure: &Structure,
    waveforms: &Waveforms,
    sampling: Sampling,
    line: &Line,
) -> Result<greens::Greens, Error> {
    let s = synthesize(structure, waveforms, sampling, line, &BASE_ELEMENTS)?;
    let stations = s
        .stations
        .into_iter()
        .map(|(name, depths)| greens::Station {
            name,
            depths: depths
                .into_iter()
                .map(|traces| traces.try_into().expect("one trace per element"))
                .collect(),
        })
        .collect();
    Ok(greens::Greens {
        title: structure.title.clone(),
        vp: s.vp,
        vs: s.vs,
        dt: sampling.dt,
        h0: line.h0,
        dk: line.dk,
        dip: line.dip,
        nk: line.nk,
        k0: line.k0,
        nt: sampling.nt,
        stations,
    })
}

/// `green3`: one mechanism, for `fort.33`.
pub fn fixed(
    structure: &Structure,
    waveforms: &Waveforms,
    sampling: Sampling,
    line: &Line,
    mechanism: Mechanism,
) -> Result<fixed::Greens, Error> {
    let s = synthesize(structure, waveforms, sampling, line, &[mechanism])?;
    let stations = s
        .stations
        .into_iter()
        .map(|(name, depths)| fixed::Station {
            name,
            depths: depths.into_iter().flatten().collect(),
        })
        .collect();
    Ok(fixed::Greens {
        title: structure.title.clone(),
        vp: s.vp,
        vs: s.vs,
        nt: sampling.nt,
        dt: sampling.dt,
        strike: mechanism.strike,
        dip: mechanism.dip,
        rake: mechanism.rake,
        h0: line.h0,
        dk: line.dk,
        nk: line.nk,
        k0: line.k0,
        stations,
    })
}

struct Synthetics {
    /// P and S velocity of the source layer at the reference depth.
    vp: f32,
    vs: f32,
    /// Each station's name and, per source point, one trace per mechanism.
    stations: Vec<(String, Vec<Vec<Vec<f32>>>)>,
}

fn synthesize(
    structure: &Structure,
    waveforms: &Waveforms,
    sampling: Sampling,
    line: &Line,
    mechanisms: &[Mechanism],
) -> Result<Synthetics, Error> {
    let nt = sampling.nt;
    let half = nt / 2;
    let df = 1.0 / (sampling.dt * nt as f32);
    let dw = (df * 2.0) * PI;
    let mut model = Model {
        source: structure.source.clone(),
        receiver: &structure.receiver,
        bounce: &structure.bounce,
        zqp: response::qfm(nt, structure.tq_p, df),
        zqs: response::qfm(nt, structure.tq_s, df),
        nt,
        dt: sampling.dt,
        df,
        dw,
    };
    // The half-space has no thickness. The source-layer search relies on it.
    if let Some(last) = model.source.last_mut() {
        last.thickness = 0.0;
    }
    let (ll, _) = source_layer(&model.source, line.h0);
    let mut out = Synthetics {
        vp: model.source[ll].vp,
        vs: model.source[ll].vs,
        stations: Vec::new(),
    };
    for station in &waveforms.stations {
        if is_skipped(station) {
            continue;
        }
        let name = short_name(station);
        let phase = Phase::from_code(station.phase).ok_or_else(|| Error {
            station: name.clone(),
            message: format!("phase code {} is not 1 to 4", station.phase),
        })?;
        let (zi, ds) = match &station.instrument {
            Instrument::Gdsn {
                power,
                magnification,
                normalisation,
                poles,
                zeros,
            } => {
                let complex = |v: &[(f32, f32)]| -> Vec<Complex> {
                    v.iter().map(|&(re, im)| Complex::new(re, im)).collect()
                };
                let zi = response::instg(
                    nt,
                    dw,
                    &complex(poles),
                    &complex(zeros),
                    *normalisation,
                    *power,
                );
                (zi, magnification * 1.0e4)
            }
            Instrument::Wwssn {
                magnification,
                seismometer_period,
                galvanometer_period,
                seismometer_damping,
                galvanometer_damping,
                coupling,
            } => {
                let zi = response::instw(
                    nt,
                    dw,
                    *seismometer_period,
                    *galvanometer_period,
                    *seismometer_damping,
                    *galvanometer_damping,
                    *coupling,
                );
                (zi, *magnification)
            }
        };
        let mut fc = station.spreading * ds;
        if phase.is_p() {
            match station.component {
                2 => fc = -(fc * (station.back_azimuth * RAD).cos()),
                3 => fc = -(fc * (station.back_azimuth * RAD).sin()),
                _ => {}
            }
        }
        let ray = Ray {
            phase,
            component: station.component,
            azimuth: station.azimuth,
            p: station.slowness,
            fc,
            zi: &zi,
        };
        let depths = (1..=line.nk)
            .map(|k| model.depth(&ray, line.depth(k), mechanisms, half))
            .collect();
        out.stations.push((name, depths));
    }
    Ok(out)
}

/// The layer holding depth `h` and the depth of its top. A depth on an
/// interface belongs to the layer above it; past every interface it is the
/// half-space.
fn source_layer(layers: &[Layer], h: f32) -> (usize, f32) {
    let mut hl = 0.0;
    for (l, layer) in layers[..layers.len() - 1].iter().enumerate() {
        hl += layer.thickness;
        if h - hl < 0.0 {
            return (l, hl - layer.thickness);
        }
    }
    let last = layers.len() - 1;
    (last, hl - layers[last].thickness)
}

/// What every station shares.
struct Model<'a> {
    source: Vec<Layer>,
    receiver: &'a [Layer],
    bounce: &'a [Layer],
    zqp: Vec<Complex>,
    zqs: Vec<Complex>,
    nt: usize,
    dt: f32,
    df: f32,
    dw: f32,
}

/// One station's ray and its gain.
struct Ray<'a> {
    phase: Phase,
    component: i32,
    azimuth: f32,
    p: f32,
    fc: f32,
    zi: &'a [Complex],
}

/// The structure's part of the synthetics at one source depth, shared by
/// the six elements.
struct Paths {
    ll: usize,
    source: Source,
    receiver: Vec<Complex>,
    /// The bounce-point reflection, for PP only.
    bounce: Vec<Complex>,
    tp: f32,
    ts: f32,
    delay: f32,
}

impl Model<'_> {
    fn depth(&self, ray: &Ray, h: f32, mechanisms: &[Mechanism], half: usize) -> Vec<Vec<f32>> {
        let paths = self.paths(ray, h, half);
        let layer = &self.source[paths.ll];
        mechanisms
            .iter()
            .map(|m| {
                let pattern = radp(
                    ray.azimuth - m.strike,
                    m.dip,
                    m.rake,
                    ray.p,
                    layer.vp,
                    layer.vs,
                );
                self.bodyw(ray, &paths, &pattern, half)
            })
            .collect()
    }

    fn paths(&self, ray: &Ray, h: f32, half: usize) -> Paths {
        let (ll, hl) = source_layer(&self.source, h);
        let dh = h - hl;
        let (p, dw) = (ray.p, self.dw);
        let (source, receiver, bounce) = if ray.phase == Phase::Sh {
            let source = propagator::reflsh(&self.source, ll, p, dw, half);
            let receiver = propagator::cnvrsh(self.receiver, p, dw, half);
            (source, receiver, Vec::new())
        } else {
            let p_wave = ray.phase.is_p();
            let source = propagator::refl(&self.source, ll, p_wave, p, dw, half);
            let receiver = propagator::cnvr(self.receiver, p_wave, p, dw, half);
            let bounce = if ray.phase == Phase::Pp {
                let last = self.bounce.len() - 1;
                propagator::refl(self.bounce, last, true, p, dw, half).pu
            } else {
                Vec::new()
            };
            (source, receiver, bounce)
        };
        let Receiver {
            horizontal,
            vertical,
            time: receiver_time,
        } = receiver;
        let receiver = if ray.phase != Phase::Sh && ray.component == 1 {
            vertical
        } else {
            horizontal
        };
        let layer = &self.source[ll];
        let p2 = p * p;
        let ts = (1.0 / (layer.vs * layer.vs) - p2).sqrt() * dh;
        // The reference computes tp for P and PP only, so its SV synthetic
        // uses whatever tp the previous call left. Here SV gets its own.
        let tp = (1.0 / (layer.vp * layer.vp) - p2).sqrt() * dh;
        let arrivals = source.time + receiver_time;
        let delay = match ray.phase {
            Phase::P => arrivals - tp,
            Phase::Pp => (arrivals - tp) - 10.0,
            Phase::Sv | Phase::Sh => arrivals - ts,
        };
        Paths {
            ll,
            source,
            receiver,
            bounce,
            tp,
            ts,
            delay,
        }
    }

    /// One synthetic trace: `sub.bodyw.f`.
    fn bodyw(&self, ray: &Ray, paths: &Paths, r: &Radiation, half: usize) -> Vec<f32> {
        let n = self.nt;
        let nl = self.source.len();
        let half_space = &self.source[nl - 1];
        let source_layer = &self.source[paths.ll];
        let tl = self.dt * n as f32;
        let (tp, ts) = (paths.tp, paths.ts);
        let s = &paths.source;
        let e = |t: f32| Complex::new(0.0, t).exp();
        // Depth phases that would arrive past the window are left out when
        // the source is in the half-space.
        let in_half_space = paths.ll == nl - 1;
        let fc = match ray.phase {
            Phase::P | Phase::Pp => {
                let v = half_space.vp;
                ray.fc / (((4.0 * PI) * half_space.density) * ((v * v) * v))
            }
            Phase::Sv => {
                let v = half_space.vs;
                ray.fc / (((4.0 * PI) * half_space.density) * ((v * v) * v))
            }
            Phase::Sh => {
                let v = half_space.vs;
                ray.fc / ((((4.0 * PI) * half_space.density) * (v * v)) * source_layer.vs)
            }
        };
        let mut z = vec![Complex::default(); n];
        for i in 0..half {
            let w = self.dw * i as f32;
            let mut zi = match ray.phase {
                Phase::P | Phase::Pp => {
                    let mut z = (s.pd[i] * c(r.pd)) * e(w * tp);
                    if !(in_half_space && tp * 2.0 >= tl) {
                        z = ((z + (s.pu[i] * c(r.pu)) * e(-(w * tp)))
                            - (s.sd[i] * c(r.svd)) * e(w * ts))
                            - (s.su[i] * c(r.svu)) * e(-(w * ts));
                    }
                    z = z * self.zqp[i];
                    if ray.phase == Phase::Pp {
                        z = ((z * paths.bounce[i]) * self.zqp[i]) * Complex::new(0.0, 1.0);
                    }
                    z
                }
                Phase::Sv => {
                    let mut z = (s.sd[i] * c(r.svd)) * e(w * ts);
                    if !(in_half_space && ts * 2.0 >= tl) {
                        z = ((z - (s.pu[i] * c(r.pu)) * e(-(w * tp)))
                            - (s.pd[i] * c(r.pd)) * e(w * tp))
                            + (s.su[i] * c(r.svu)) * e(-(w * ts));
                    }
                    z * self.zqs[i]
                }
                Phase::Sh => {
                    let mut z = (s.sd[i] * c(r.shd)) * e(w * ts);
                    if !(in_half_space && ts * 2.0 >= tl) {
                        z = z + (s.su[i] * c(r.shu)) * e(-(w * ts));
                    }
                    z * self.zqs[i]
                }
            };
            zi = (((c(fc) * zi) * paths.receiver[i]) * ray.zi[i]) * e(w * paths.delay);
            z[i] = zi;
            if i > 0 {
                z[n - i] = zi.conj();
            }
        }
        z[half] = Complex::default();
        cfft(&mut z, true);
        // The full product keeps the sign of a zero as the reference does.
        z.iter().map(|v| (*v * c(self.df)).re).collect()
    }
}

/// A real operand in a complex expression.
fn c(re: f32) -> Complex {
    Complex::from(re)
}
