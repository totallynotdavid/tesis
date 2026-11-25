//! `fort.33`: the Green's functions of one fixed mechanism, which `green3`
//! writes and `inversion3` reads, as a gfortran unformatted file.
//!
//! The header record holds the structure title, the source layer velocities
//! at the reference depth, the sampling, the mechanism and the depth grid.
//! Each station then has its name record and one record per sample, depth
//! outer.

use super::greens::{NAME, TITLE};
use super::unformatted::{Error, Reader, padded, reals, record, word};

#[derive(Debug, Clone, PartialEq)]
pub struct Greens {
    /// The `fort.2` title, padded or cut to 40 characters.
    pub title: String,
    /// P and S velocity of the source layer at the reference depth `h0`.
    pub vp: f32,
    pub vs: f32,
    /// Samples per trace.
    pub nt: usize,
    pub dt: f32,
    /// The mechanism in degrees.
    pub strike: f32,
    pub dip: f32,
    pub rake: f32,
    /// The reference depth and the spacing of the grid points along the
    /// fault dip, in km.
    pub h0: f32,
    pub dk: f32,
    /// Number of grid points, and the one at `h0`, counted from 1.
    pub nk: usize,
    pub k0: i32,
    pub stations: Vec<Station>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Station {
    /// The first 12 characters of the `fort.1` name record, padded.
    pub name: String,
    /// One trace per grid point.
    pub depths: Vec<Vec<f32>>,
}

pub fn write(greens: &Greens) -> Vec<u8> {
    let mut out = Vec::new();
    let mut header = padded(&greens.title, TITLE);
    header.extend(reals(&[greens.vp, greens.vs]));
    header.extend((greens.nt as i32).to_le_bytes());
    header.extend(reals(&[
        greens.dt,
        greens.strike,
        greens.dip,
        greens.rake,
        greens.h0,
        greens.dk,
    ]));
    for v in [greens.nk as i32, greens.k0] {
        header.extend(v.to_le_bytes());
    }
    record(&mut out, &header);
    for station in &greens.stations {
        record(&mut out, &padded(&station.name, NAME));
        for sample in station.depths.iter().flatten() {
            record(&mut out, &sample.to_le_bytes());
        }
    }
    out
}

pub fn read(bytes: &[u8]) -> Result<Greens, Error> {
    let mut records = Reader::new(bytes);
    let header = records.next(TITLE + 11 * 4)?;
    let real = |i: usize| f32::from_le_bytes(word(&header[TITLE + 4 * i..]));
    let int = |i: usize| i32::from_le_bytes(word(&header[TITLE + 4 * i..]));
    let (nt, nk) = (records.count(int(2))?, records.count(int(9))?);
    let mut greens = Greens {
        title: String::from_utf8_lossy(&header[..TITLE]).into_owned(),
        vp: real(0),
        vs: real(1),
        nt,
        dt: real(3),
        strike: real(4),
        dip: real(5),
        rake: real(6),
        h0: real(7),
        dk: real(8),
        nk,
        k0: int(10),
        stations: Vec::new(),
    };
    while !records.at_end() {
        let name = records.next_text(NAME)?;
        let mut depths = Vec::with_capacity(nk);
        for _ in 0..nk {
            let mut trace = Vec::with_capacity(nt);
            for _ in 0..nt {
                trace.extend(records.next_reals(1)?);
            }
            depths.push(trace);
        }
        greens.stations.push(Station { name, depths });
    }
    Ok(greens)
}
