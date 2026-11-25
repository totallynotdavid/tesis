//! `fort.1`: the event title, then one block per station with its geometry,
//! instrument and normalised waveform. The field meanings follow section 1 of
//! `reference/manual.pdf`.

use std::fmt::Write as _;

use super::records::{Error, Records};
use super::{reals, write_rows};

#[derive(Debug, Clone, PartialEq)]
pub struct Waveforms {
    pub title: String,
    pub stations: Vec<Station>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Station {
    /// The whole name record. The programs read its first 12 characters;
    /// `conv.sac.farm` writes picking notes after them.
    pub name: String,
    pub azimuth: f32,
    pub back_azimuth: f32,
    /// Epicentral distance in degrees.
    pub distance: f32,
    /// Ray parameter in s/km.
    pub slowness: f32,
    /// Geometrical spreading factor.
    pub spreading: f32,
    /// Samples before the P or S onset (`ix0`). A station with `ix0` past
    /// its last sample is skipped by `green`.
    pub onset: i32,
    /// `ib`: 1 P, 2 SV, 3 SH, 4 PP.
    pub phase: i32,
    /// `ic`: 1 Z, 2 N, 3 E. Ignored for SH; for SV, 1 vertical, 2 horizontal.
    pub component: i32,
    pub instrument: Instrument,
    /// Maximum amplitude, in cm for WWSSN and micrometres for GDSN.
    pub peak: f32,
    /// Sample interval in seconds.
    pub dt: f32,
    /// Samples divided by `peak`.
    pub samples: Vec<f32>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Instrument {
    /// `im = 0`: a WWSSN long-period response.
    Wwssn {
        magnification: f32,
        seismometer_period: f32,
        galvanometer_period: f32,
        seismometer_damping: f32,
        galvanometer_damping: f32,
        coupling: f32,
    },
    /// `im = 1`: a response in poles and zeros. With none, `power` 0 means
    /// displacement and 1 velocity.
    Gdsn {
        power: i32,
        magnification: f32,
        normalisation: f32,
        poles: Vec<(f32, f32)>,
        zeros: Vec<(f32, f32)>,
    },
}

pub fn parse(text: &str) -> Result<Waveforms, Error> {
    let mut records = Records::new(text);
    let title = records.raw().unwrap_or_default().trim_end().to_owned();
    let mut stations = Vec::new();
    while let Some(name) = records.raw() {
        if name.trim().is_empty() {
            continue;
        }
        stations.push(parse_station(name, &mut records)?);
    }
    Ok(Waveforms { title, stations })
}

fn parse_station(name: &str, records: &mut Records) -> Result<Station, Error> {
    let mut read = records.read();
    let (azimuth, back_azimuth, distance) = (read.f32()?, read.f32()?, read.f32()?);
    let (slowness, spreading, onset) = (read.f32()?, read.f32()?, read.i32()?);
    drop(read);
    let mut read = records.read();
    let (im, phase, component) = (read.i32()?, read.i32()?, read.i32()?);
    drop(read);
    let instrument = if im == 1 {
        let mut read = records.read();
        let (n_poles, n_zeros, power) = (read.count()?, read.count()?, read.i32()?);
        drop(read);
        let mut read = records.read();
        let (magnification, normalisation) = (read.f32()?, read.f32()?);
        drop(read);
        // One pole or zero per record: `inversion` skips them a record each.
        let mut pairs = |n: usize| -> Result<Vec<(f32, f32)>, Error> {
            (0..n)
                .map(|_| {
                    let mut read = records.read();
                    Ok((read.f32()?, read.f32()?))
                })
                .collect()
        };
        let poles = pairs(n_poles)?;
        let zeros = pairs(n_zeros)?;
        Instrument::Gdsn {
            power,
            magnification,
            normalisation,
            poles,
            zeros,
        }
    } else {
        let v = records.read().f32s(6)?;
        Instrument::Wwssn {
            magnification: v[0],
            seismometer_period: v[1],
            galvanometer_period: v[2],
            seismometer_damping: v[3],
            galvanometer_damping: v[4],
            coupling: v[5],
        }
    };
    let mut read = records.read();
    let (peak, count, dt) = (read.f32()?, read.count()?, read.f32()?);
    drop(read);
    let samples = records.read().f32s(count)?;
    Ok(Station {
        name: name.trim_end().to_owned(),
        azimuth,
        back_azimuth,
        distance,
        slowness,
        spreading,
        onset,
        phase,
        component,
        instrument,
        peak,
        dt,
        samples,
    })
}

/// Writes `fort.1`. Each real is the shortest decimal that reads back to the
/// same `f32`, so a Fortran reader gets the values bit for bit.
pub fn write(waveforms: &Waveforms) -> String {
    let mut out = format!("{}\n", waveforms.title);
    for s in &waveforms.stations {
        let _ = writeln!(out, "{}", s.name);
        let geometry = [
            s.azimuth,
            s.back_azimuth,
            s.distance,
            s.slowness,
            s.spreading,
        ];
        let _ = writeln!(out, "{} {:>5}", reals(&geometry), s.onset);
        match &s.instrument {
            Instrument::Gdsn {
                power,
                magnification,
                normalisation,
                poles,
                zeros,
            } => {
                let _ = writeln!(out, "  1 {:>2} {:>2}", s.phase, s.component);
                let _ = writeln!(out, " {:>4} {:>4} {:>4}", poles.len(), zeros.len(), power);
                let _ = writeln!(out, "{}", reals(&[*magnification, *normalisation]));
                for (re, im) in poles.iter().chain(zeros) {
                    let _ = writeln!(out, "{}", reals(&[*re, *im]));
                }
            }
            Instrument::Wwssn {
                magnification,
                seismometer_period,
                galvanometer_period,
                seismometer_damping,
                galvanometer_damping,
                coupling,
            } => {
                let _ = writeln!(out, "  0 {:>2} {:>2}", s.phase, s.component);
                let values = [
                    *magnification,
                    *seismometer_period,
                    *galvanometer_period,
                    *seismometer_damping,
                    *galvanometer_damping,
                    *coupling,
                ];
                let _ = writeln!(out, "{}", reals(&values));
            }
        }
        let _ = writeln!(out, " {:>11} {:>4} {:>7}", s.peak, s.samples.len(), s.dt);
        write_rows(&mut out, &s.samples, 10);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const GDSN: &str = "Test event
ABC         notes after column 12
  10.5 -170.0 45.0 0.07 12.8 3
  1 1 1
  2 1 0
  1.0 1.0
  -0.1 0.2
  -0.1 -0.2
  0.0 0.0
  0.5e+02 3 0.5
 0.1 -0.2
 1.0
";

    #[test]
    fn a_gdsn_block_reads_its_poles_and_zeros() {
        let station = &parse(GDSN).unwrap().stations[0];
        assert_eq!(station.name, "ABC         notes after column 12");
        assert_eq!(station.onset, 3);
        let Instrument::Gdsn { poles, zeros, .. } = &station.instrument else {
            panic!("expected GDSN");
        };
        assert_eq!(poles, &[(-0.1, 0.2), (-0.1, -0.2)]);
        assert_eq!(zeros, &[(0.0, 0.0)]);
        assert_eq!(station.samples, [0.1, -0.2, 1.0]);
    }

    #[test]
    fn writing_then_reading_gives_the_same_values() {
        let waveforms = parse(GDSN).unwrap();
        assert_eq!(parse(&write(&waveforms)).unwrap(), waveforms);
    }
}
