//! The event file: the source, which stations to fetch, and how to cut and
//! filter their P windows. The fetched files live beside it.

use std::path::{Path, PathBuf};

use libmseed::NsTime;
use serde::{Deserialize, Deserializer};

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
}

impl Files {
    /// The files beside the event file at `path`.
    pub fn beside(path: &Path) -> Self {
        let dir = path.parent().unwrap_or(Path::new("."));
        Self {
            stations: dir.join("stations.xml"),
            waveforms: dir.join("waveforms.mseed"),
            fort1: dir.join("fort.1"),
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
