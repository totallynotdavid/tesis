//! `fort.3`: the Green's functions of the six base elements, which `green`
//! writes and `inversion` reads, as a gfortran unformatted file.
//!
//! The header record holds the structure title, the source layer velocities
//! at the reference depth, the sampling and the depth grid. Each station then
//! has its name record and one record of `nt` samples per depth and base
//! element, depth outer.

use super::unformatted::{Error, Reader, padded, reals, record, word};

/// The six base elements per depth, in `green`'s order: Mxy, Mzz-Myy, Myz,
/// Mzx, Mzz-Mxx and the isotropic Mxx+Myy+Mzz.
pub const ELEMENTS: usize = 6;

pub(super) const TITLE: usize = 40;
pub(super) const NAME: usize = 12;

#[derive(Debug, Clone, PartialEq)]
pub struct Greens {
    /// The `fort.2` title, padded or cut to 40 characters.
    pub title: String,
    /// P and S velocity of the source layer at the reference depth `h0`.
    pub vp: f32,
    pub vs: f32,
    pub dt: f32,
    /// The reference depth and the spacing of the depth grid, in km.
    pub h0: f32,
    pub dk: f32,
    /// Dip of the line of grid points, in degrees.
    pub dip: f32,
    /// Number of grid points, and the one at `h0`, counted from 1.
    pub nk: usize,
    pub k0: i32,
    /// Samples per trace.
    pub nt: usize,
    pub stations: Vec<Station>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Station {
    /// The first 12 characters of the `fort.1` name record, padded.
    pub name: String,
    /// One trace per base element for each grid point.
    pub depths: Vec<[Vec<f32>; ELEMENTS]>,
}

pub fn write(greens: &Greens) -> Vec<u8> {
    let mut out = Vec::new();
    let mut header = padded(&greens.title, TITLE);
    header.extend(reals(&[
        greens.vp, greens.vs, greens.dt, greens.h0, greens.dk, greens.dip,
    ]));
    for v in [greens.nk as i32, greens.k0, greens.nt as i32] {
        header.extend(v.to_le_bytes());
    }
    record(&mut out, &header);
    for station in &greens.stations {
        record(&mut out, &padded(&station.name, NAME));
        for trace in station.depths.iter().flatten() {
            record(&mut out, &reals(trace));
        }
    }
    out
}

pub fn read(bytes: &[u8]) -> Result<Greens, Error> {
    let mut records = Reader::new(bytes);
    let header = records.next(TITLE + 9 * 4)?;
    let real = |i: usize| f32::from_le_bytes(word(&header[TITLE + 4 * i..]));
    let int = |i: usize| i32::from_le_bytes(word(&header[TITLE + 4 * i..]));
    let (nk, nt) = (records.count(int(6))?, records.count(int(8))?);
    let mut greens = Greens {
        title: String::from_utf8_lossy(&header[..TITLE]).into_owned(),
        vp: real(0),
        vs: real(1),
        dt: real(2),
        h0: real(3),
        dk: real(4),
        dip: real(5),
        nk,
        k0: int(7),
        nt,
        stations: Vec::new(),
    };
    while !records.at_end() {
        let name = records.next_text(NAME)?;
        let mut depths = Vec::with_capacity(nk);
        for _ in 0..nk {
            let mut traces: [Vec<f32>; ELEMENTS] = Default::default();
            for trace in &mut traces {
                *trace = records.next_reals(nt)?;
            }
            depths.push(traces);
        }
        greens.stations.push(Station { name, depths });
    }
    Ok(greens)
}
