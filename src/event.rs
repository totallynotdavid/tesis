//! The event file: the source, which stations to fetch, and how to cut and
//! filter their P windows. The fetched files live beside it.

use std::path::{Path, PathBuf};

use libmseed::NsTime;
use serde::{Deserialize, Deserializer};

use crate::formats::structure::{Layer, Structure};
use crate::green::{Line, Sampling};
use crate::prep::earth::Phase;
use crate::prep::{Event, Quantity, Time, Window};

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    /// The `fort.1` title; the programs keep its first 20 characters.
    pub title: String,
    /// UTC, ISO 8601.
    #[serde(deserialize_with = "time")]
    pub origin: NsTime,
    pub latitude: f32,
    pub longitude: f32,
    /// Kilometres.
    pub depth: f32,
    #[serde(default)]
    pub output: Quantity,
    /// Degrees.
    #[serde(default = "distances")]
    pub distances: [f32; 2],
    #[serde(default)]
    pub stations: Stations,
    pub p: P,
    pub structure: Option<Model>,
    #[serde(default)]
    pub green: Green,
}

/// The `fort.2` velocity structure. Each layer is `[vp, vs, density,
/// thickness]` in km/s, g/cm3 and km, top down; the last is the half-space.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Model {
    pub title: String,
    pub tq_p: f32,
    pub tq_s: f32,
    pub source: Vec<[f32; 4]>,
    pub receiver: Vec<[f32; 4]>,
    pub bounce: Vec<[f32; 4]>,
}

/// The sampling of the Green's functions and the line of source points.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct Green {
    pub nt: usize,
    /// Defaults to `p.dt`.
    pub dt: Option<f32>,
    /// Depth of the reference point in km; defaults to the event depth.
    pub h0: Option<f32>,
    pub spacing: f32,
    pub points: usize,
    /// The point at `h0`, counted from 1.
    pub reference: i32,
    pub dip: f32,
}

impl Default for Green {
    fn default() -> Self {
        Self {
            nt: 256,
            dt: None,
            h0: None,
            spacing: 5.0,
            points: 1,
            reference: 1,
            dip: 90.0,
        }
    }
}

/// An FDSN web service and the channels to ask it for, in its query syntax.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct Stations {
    pub service: String,
    pub networks: String,
    pub location: String,
    pub channel: String,
}

impl Default for Stations {
    fn default() -> Self {
        Self {
            service: "https://service.earthscope.org".into(),
            networks: "IU,II".into(),
            location: "00".into(),
            channel: "BHZ".into(),
        }
    }
}

/// The P window. The fields mean what [`Window`]'s do.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct P {
    pub lead: f32,
    pub duration: f32,
    pub dt: f32,
    #[serde(default)]
    pub low: f32,
    #[serde(default)]
    pub high: f32,
}

/// The paths of an event's files.
pub struct Files {
    pub stations: PathBuf,
    pub waveforms: PathBuf,
    pub fort1: PathBuf,
    pub fort3: PathBuf,
}

impl Config {
    pub fn parse(text: &str) -> Result<Self, String> {
        toml::from_str(text).map_err(|e| e.to_string())
    }

    pub fn event(&self) -> Event {
        Event {
            title: self.title.clone(),
            latitude: self.latitude,
            longitude: self.longitude,
            origin: Time::from_nanos(self.origin.as_nanos()),
            output: self.output,
            distances: self.distances[0]..=self.distances[1],
        }
    }

    pub fn window(&self) -> Window {
        Window {
            phase: Phase::P,
            lead: self.p.lead,
            duration: self.p.duration,
            dt: self.p.dt,
            low: self.p.low,
            high: self.p.high,
        }
    }

    pub fn structure(&self) -> Result<Structure, String> {
        let model = self.structure.as_ref().ok_or("no [structure] table")?;
        let layers = |rows: &[[f32; 4]], name: &str| -> Result<Vec<Layer>, String> {
            if rows.is_empty() {
                return Err(format!("structure.{name} has no layers"));
            }
            Ok(rows
                .iter()
                .map(|&[vp, vs, density, thickness]| Layer {
                    vp,
                    vs,
                    density,
                    thickness,
                })
                .collect())
        };
        Ok(Structure {
            title: model.title.clone(),
            tq_p: model.tq_p,
            tq_s: model.tq_s,
            source: layers(&model.source, "source")?,
            receiver: layers(&model.receiver, "receiver")?,
            bounce: layers(&model.bounce, "bounce")?,
        })
    }

    pub fn sampling(&self) -> Result<Sampling, String> {
        let nt = self.green.nt;
        if nt < 2 || !nt.is_power_of_two() {
            return Err(format!("green.nt is {nt}, not a power of two"));
        }
        Ok(Sampling {
            nt,
            dt: self.green.dt.unwrap_or(self.p.dt),
        })
    }

    pub fn line(&self) -> Line {
        Line {
            h0: self.green.h0.unwrap_or(self.depth),
            dk: self.green.spacing,
            nk: self.green.points,
            k0: self.green.reference,
            dip: self.green.dip,
        }
    }
}

impl Files {
    /// The files beside the event file at `path`.
    pub fn beside(path: &Path) -> Self {
        let dir = path.parent().unwrap_or(Path::new("."));
        Self {
            stations: dir.join("stations.xml"),
            waveforms: dir.join("waveforms.mseed"),
            fort1: dir.join("fort.1"),
            fort3: dir.join("fort.3"),
        }
    }
}

fn distances() -> [f32; 2] {
    [30.0, 90.0]
}

fn time<'de, D: Deserializer<'de>>(d: D) -> Result<NsTime, D::Error> {
    let text = String::deserialize(d)?;
    NsTime::parse(&text).map_err(serde::de::Error::custom)
}
