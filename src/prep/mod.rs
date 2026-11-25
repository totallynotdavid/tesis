//! Turns a raw seismogram into a `fort.1` station block, as `conv.sac.farm`
//! does: cut a window around the J-B onset, remove the instrument response
//! and band-pass in the frequency domain, decimate, and normalise.

pub mod earth;

use std::fmt;
use std::ops::RangeInclusive;

use crate::fft::{Complex, cfft, imag_minus};
use crate::formats::waveforms::{Instrument, Station};
use earth::{Earth, Phase};

#[expect(
    clippy::approx_constant,
    reason = "the reference's value, one ulp below f32::consts::PI"
)]
const PI: f32 = 3.1415926;

/// A ground-motion quantity. The number is how many times it differentiates
/// displacement.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Quantity {
    #[default]
    Displacement = 0,
    Velocity = 1,
    Acceleration = 2,
}

/// A time as minutes since some epoch and seconds within the minute. The
/// epoch only has to be common to the event and its records.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Time {
    pub minute: i64,
    pub second: f32,
}

impl Time {
    pub fn from_nanos(nanos: i64) -> Self {
        const MINUTE: i64 = 60_000_000_000;
        Self {
            minute: nanos.div_euclid(MINUTE),
            second: (nanos.rem_euclid(MINUTE) as f64 / 1e9) as f32,
        }
    }
}

pub struct Event {
    pub title: String,
    pub latitude: f32,
    pub longitude: f32,
    pub origin: Time,
    /// The quantity the stations are converted to.
    pub output: Quantity,
    /// Stations outside this range of distance, in degrees, are skipped.
    pub distances: RangeInclusive<f32>,
}

/// One channel's samples, from its first sample.
pub struct Record {
    pub code: String,
    /// The SEED channel code, such as `BHZ`.
    pub channel: String,
    pub latitude: f32,
    pub longitude: f32,
    pub start: Time,
    pub dt: f32,
    pub samples: Vec<f32>,
    pub response: Response,
}

/// Poles and zeros in rad/s, from ground motion in metres to counts.
#[derive(Debug, Clone)]
pub struct Response {
    /// The ground motion the poles and zeros take as input.
    pub units: Quantity,
    pub normalisation: f32,
    /// Counts per unit of ground motion.
    pub sensitivity: f32,
    pub poles: Vec<Complex>,
    pub zeros: Vec<Complex>,
}

/// What to cut and how to filter it.
pub struct Window {
    pub phase: Phase,
    /// Seconds of record kept before the onset.
    pub lead: f32,
    /// Seconds of record kept from the window start.
    pub duration: f32,
    /// The output sample interval. It must be a whole multiple of the record's.
    pub dt: f32,
    /// Band-pass corners in Hz; 0 leaves that side open.
    pub low: f32,
    pub high: f32,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Skip {
    Distance(f32),
    Channel(String),
    /// The window starts before the first sample.
    BeforeRecord,
    /// The window starts after the last sample.
    AfterRecord,
    ShortWindow,
    /// The output interval is shorter than the record's.
    Upsampling,
    /// The window is zero after filtering.
    Flat,
}

impl fmt::Display for Skip {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Distance(d) => write!(f, "distance {d:.1} degrees is out of range"),
            Self::Channel(c) => write!(f, "channel {c:?} has no Z, N, E, 1 or 2 component"),
            Self::BeforeRecord => f.write_str("the window starts before the record"),
            Self::AfterRecord => f.write_str("the window starts after the record"),
            Self::ShortWindow => f.write_str("the window holds no samples"),
            Self::Upsampling => f.write_str("the output interval is shorter than the record's"),
            Self::Flat => f.write_str("the filtered window is zero"),
        }
    }
}

/// The `fort.1` title line.
pub fn title(event: &Event, earth: &Earth) -> String {
    let quantity = ["disp.", "vel. ", "acc. "][event.output as usize];
    let hour = event.origin.minute.div_euclid(60).rem_euclid(24);
    let minute = event.origin.minute.rem_euclid(60);
    format!(
        "{:<20.20}{quantity}{:>7.2}{:>7.2}{:>7.2}{hour:>3}{minute:>3}{:>6.2}",
        event.title,
        event.latitude,
        event.longitude,
        earth.depth(),
        event.origin.second,
    )
}

/// The `fort.1` block for one record.
pub fn convert(
    event: &Event,
    earth: &Earth,
    record: &Record,
    window: &Window,
) -> Result<Station, Skip> {
    let component =
        component(&record.channel).ok_or_else(|| Skip::Channel(record.channel.clone()))?;
    let geometry = earth.geometry(
        (event.latitude, event.longitude),
        (record.latitude, record.longitude),
    );
    if !event.distances.contains(&geometry.distance) {
        return Err(Skip::Distance(geometry.distance));
    }
    let travel = earth.travel(window.phase, geometry.distance);
    let y = cut(event, record, window, travel)?;
    let y = deconvolve(&y, record, window, event.output);
    let peak = y.iter().fold(0.0f32, |m, v| v.abs().max(m));
    if peak == 0.0 {
        return Err(Skip::Flat);
    }
    let step = (window.dt / record.dt + 0.1) as usize;
    if step < 1 {
        return Err(Skip::Upsampling);
    }
    let samples = (0..y.len() / step).map(|i| y[i * step] / peak).collect();
    Ok(Station {
        name: name(event, record, window, travel),
        azimuth: geometry.azimuth,
        back_azimuth: geometry.back_azimuth,
        distance: geometry.distance,
        slowness: if geometry.distance < 10.0 {
            0.0
        } else {
            geometry.slowness
        },
        spreading: geometry.spreading,
        onset: (window.lead / window.dt + 0.5) as i32,
        phase: match window.phase {
            Phase::P => 1,
            Phase::S => 2,
        },
        component,
        instrument: Instrument::Gdsn {
            power: event.output as i32,
            magnification: 1.0,
            normalisation: 1.0,
            poles: Vec::new(),
            zeros: Vec::new(),
        },
        peak,
        dt: window.dt,
        samples,
    })
}

/// `ic`: 1 vertical, 2 north, 3 east. Channels ending in 1 and 2 count as
/// east and north.
fn component(channel: &str) -> Option<i32> {
    match channel.chars().last()? {
        'Z' => Some(1),
        'N' | '2' => Some(2),
        'E' | '1' => Some(3),
        _ => None,
    }
}

/// The window's samples with their mean removed. `travel` is the onset in
/// seconds after the origin.
fn cut(event: &Event, record: &Record, window: &Window, travel: f32) -> Result<Vec<f32>, Skip> {
    let minutes = (60 * (event.origin.minute - record.start.minute)) as f32;
    let start = minutes + event.origin.second + travel - record.start.second - window.lead;
    if start < 0.0 {
        return Err(Skip::BeforeRecord);
    }
    let n = record.samples.len();
    let first = (start / record.dt + 1.5) as usize;
    if first > n {
        return Err(Skip::AfterRecord);
    }
    let last = (((start + window.duration) / record.dt + 1.5) as usize).min(n);
    if last < first {
        return Err(Skip::ShortWindow);
    }
    let y = &record.samples[first - 1..last];
    let mut sum = 0.0f32;
    for v in y {
        sum += v;
    }
    let mean = sum / y.len() as f32;
    Ok(y.iter().map(|v| v - mean).collect())
}

/// Divides the instrument response out of `y`, band-passes it with zero
/// phase, and scales it to micrometres (per second, per second squared) of
/// `output`.
fn deconvolve(y: &[f32], record: &Record, window: &Window, output: Quantity) -> Vec<f32> {
    let response = &record.response;
    let nn = y.len().next_power_of_two();
    let mut z = vec![Complex::default(); nn];
    for (z, &y) in z.iter_mut().zip(y) {
        *z = Complex::new(y, 0.0);
    }
    cfft(&mut z, false);
    let df = 1.0 / (record.dt * nn as f32);
    let dw = df * PI * 2.0;
    // Each extra zero at the origin integrates once; each factor `i 2 pi f`
    // differentiates once.
    let integrations = response.units as i32 - output as i32;
    let mut zeros = response.zeros.clone();
    zeros.resize(
        zeros.len() + integrations.max(0) as usize,
        Complex::default(),
    );
    for i in 2..=nn / 2 + 1 {
        let f = df * (i - 1) as f32;
        let w = dw * (i - 1) as f32;
        let low = if window.low == 0.0 {
            1.0
        } else {
            1.0 - 1.0 / (1.0 + power(f / window.low))
        };
        let high = if window.high == 0.0 {
            1.0
        } else {
            1.0 / (1.0 + power(f / window.high))
        };
        let mut zt = Complex::new(1.0, 0.0);
        for &p in &response.poles {
            zt = zt * imag_minus(w, p);
        }
        for &q in &zeros {
            zt = zt / imag_minus(w, q);
        }
        let mut zc = zt.scale(low * high);
        for _ in integrations..0 {
            zc = zc.mul_imag(2.0 * PI * f);
        }
        z[i - 1] = z[i - 1] * zc;
        if i == nn / 2 + 1 {
            z[i - 1].im = 0.0;
        }
        z[nn - i + 1] = z[i - 1].conj();
    }
    z[0] = Complex::default();
    cfft(&mut z, true);
    let scale = 1.0e6 / (response.normalisation * response.sensitivity);
    z[..y.len()]
        .iter()
        .map(|z| z.re / nn as f32 * scale)
        .collect()
}

/// `x^8`, the band-pass order, by repeated squaring as gfortran expands
/// `X**8`.
fn power(x: f32) -> f32 {
    let x2 = x * x;
    let x4 = x2 * x2;
    x4 * x4
}

/// The name record: the station code, then the onset as minute of the hour
/// and second, the lead, and the band-pass corners.
fn name(event: &Event, record: &Record, window: &Window, travel: f32) -> String {
    let mut second = event.origin.second + travel;
    let minutes = (second / 60.0) as i64;
    second -= (minutes * 60) as f32;
    let minute = event.origin.minute.rem_euclid(60) + minutes;
    format!(
        "{:<12.12}{minute:>4}{second:>7.2}{:>8.2}{:>7.2}{:>7.2}",
        record.code, -window.lead, window.low, window.high
    )
}
