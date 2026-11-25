//! Travel times, distances and ray factors from the Jeffreys-Bullen tables,
//! as `conv.sac.farm` computes them in `jbtravel`, `GEOM` and `parm`.
//!
//! `jb.ptime` and `jb.stime` give P and S times by distance (rows, in whole
//! degrees) and source depth (14 columns at [`DEPTHS`]). `jb.table` gives a
//! velocity profile; `GEOM` traces rays through it to get the slowness and
//! the geometrical spreading for each whole degree of distance.

use std::fmt;

use crate::formats::records::{self, Records};

/// Degrees to radians, as the reference rounds it.
const RAD: f32 = 0.017453292;

#[expect(
    clippy::approx_constant,
    reason = "the reference's value, one ulp above f32::consts::PI"
)]
const PI: f32 = 3.141593;

/// The source depths of the travel-time columns, in km.
pub const DEPTHS: [f32; 14] = [
    0.0, 33.0, 97.0, 160.0, 224.0, 288.0, 352.0, 415.0, 479.0, 543.0, 606.0, 670.0, 734.0, 798.0,
];

/// The reference sizes its tables for 500 degrees and reads unset entries as
/// zero.
const SPAN: usize = 501;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase {
    P,
    S,
}

#[derive(Debug, PartialEq)]
pub enum Error {
    Table {
        file: &'static str,
        error: records::Error,
    },
    Profile(&'static str),
    /// A source depth outside the travel-time columns.
    Depth(f32),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Table { file, error } => write!(f, "{file}: {error}"),
            Self::Profile(message) => write!(f, "jb.table: {message}"),
            Self::Depth(depth) => write!(
                f,
                "depth {depth} km is outside the J-B tables (0 to {} km)",
                DEPTHS[13]
            ),
        }
    }
}

impl std::error::Error for Error {}

pub struct Tables {
    p: Travel,
    s: Travel,
    profile: Profile,
}

/// One travel-time file: per whole degree, (minutes, seconds) per column.
struct Travel {
    rows: Vec<[(i32, f32); 14]>,
}

/// `jb.table`: velocities every `dh1` km, interpolated every `dh` km.
struct Profile {
    dh1: f64,
    dh: f64,
    nddh: usize,
    velocities: Vec<f64>,
}

impl Tables {
    /// The tables the reference distribution ships, built into the binary.
    pub fn jb() -> Self {
        Self::read(
            include_str!("../../reference/examples/input/jb.ptime"),
            include_str!("../../reference/examples/input/jb.stime"),
            include_str!("../../reference/examples/input/jb.table"),
        )
        .expect("the shipped J-B tables parse")
    }

    /// Reads the texts of `jb.ptime`, `jb.stime` and `jb.table`.
    pub fn read(ptime: &str, stime: &str, table: &str) -> Result<Self, Error> {
        let table_error = |file| move |error| Error::Table { file, error };
        let p = Travel::read(ptime).map_err(table_error("jb.ptime"))?;
        let s = Travel::read(stime).map_err(table_error("jb.stime"))?;
        let profile = Profile::read(table).map_err(table_error("jb.table"))?;
        if profile.velocities.len() < 2 || profile.nddh == 0 || profile.dh <= 0.0 {
            return Err(Error::Profile("needs two velocities and a positive step"));
        }
        Ok(Self { p, s, profile })
    }

    /// The tables for a source at `depth` km.
    pub fn at(&self, depth: f32) -> Result<Earth<'_>, Error> {
        if !(0.0..DEPTHS[13]).contains(&depth) {
            return Err(Error::Depth(depth));
        }
        let (spreading, slowness) = self.profile.rays(f64::from(depth));
        Ok(Earth {
            tables: self,
            depth,
            spreading,
            slowness,
        })
    }
}

impl Travel {
    fn read(text: &str) -> Result<Self, records::Error> {
        let mut records = Records::new(text);
        let mut rows = vec![[(0, 0.0); 14]; 105];
        loop {
            let mut read = records.read();
            if read.at_end() {
                break;
            }
            let n = read.count()?;
            let mut row = [(0, 0.0); 14];
            for cell in &mut row {
                *cell = (read.i32()?, read.f32()?);
            }
            if rows.len() <= n {
                rows.resize(n + 1, [(0, 0.0); 14]);
            }
            rows[n] = row;
        }
        Ok(Self { rows })
    }

    /// The time at whole degree `row` for a depth between columns `c` and
    /// `c + 1`.
    fn at(&self, row: usize, c: usize, ratio: f32) -> f32 {
        let r = &self.rows[row];
        let ((m0, s0), (m1, s1)) = (r[c], r[c + 1]);
        ((m0 * 60) as f32 + s0) + ratio * ((((m1 - m0) * 60) as f32 + s1) - s0)
    }
}

impl Profile {
    fn read(text: &str) -> Result<Self, records::Error> {
        let mut records = Records::new(text);
        records.raw();
        let mut read = records.read();
        let (dh1, count, dh, nddh) = (read.f64()?, read.count()?, read.f64()?, read.count()?);
        drop(read);
        let mut read = records.read();
        let velocities = (0..count).map(|_| read.f64()).collect::<Result<_, _>>()?;
        Ok(Self {
            dh1,
            dh,
            nddh,
            velocities,
        })
    }

    /// `GEOM`: the spreading factor and the slowness for each whole degree
    /// from 20 to 115, for a source at `h` km.
    fn rays(&self, h: f64) -> (Vec<f64>, Vec<f64>) {
        let (vv, nddh, dh) = (&self.velocities, self.nddh, self.dh);
        let ndep = nddh * (vv.len() - 1);
        // 1-based, as `V(I)` and `R(I)` are: the index is the layer.
        let mut v = vec![0.0; ndep + 1];
        for i in 2..=vv.len() {
            for j in 1..=nddh {
                v[nddh * (i - 2) + j] =
                    vv[i - 2] + (vv[i - 1] - vv[i - 2]) / self.dh1 * dh * j as f64;
            }
        }
        let r: Vec<f32> = (0..=ndep)
            .map(|i| (6371.0 - dh * i as f64) as f32)
            .collect();
        let l = (h / dh + 1.0) as usize;
        let rh = (6371.0 - h) as f32;

        let mut delta = vec![0.0; SPAN];
        for (j, slot) in delta.iter_mut().enumerate().take(181).skip(1) {
            let th = f64::from(RAD * j as f32 * 0.5);
            let p = f64::from(rh) * th.sin() / v[l];
            *slot = if p * v[ndep] / f64::from(r[ndep]) <= 1.0 {
                f64::from(PI)
            } else {
                ray_distance(p, l, rh, &v, &r)
            };
        }

        let mut tih = vec![0.0; 131];
        for (j, slot) in tih.iter_mut().enumerate().take(126).skip(5) {
            let del = f64::from(RAD * j as f32);
            let mut i = 1;
            while del <= delta[i] {
                i += 1;
            }
            let m = i - 1;
            *slot = (m as f64 + (delta[m] - del) / (delta[m] - delta[m + 1])) * 0.5;
        }

        let mut spreading = vec![0.0; SPAN];
        let mut slowness = vec![0.0; SPAN];
        for j in 20..=115 {
            let mut dtih = 0.0;
            for kk in -15i32..=15 {
                dtih += f64::from(kk) * tih[(j as i32 - kk) as usize];
            }
            dtih /= 2480.0;
            let sh = (tih[j] * f64::from(RAD)).sin();
            slowness[j] = sh / v[l];
            let pp = v[1] * slowness[j];
            let cs0 = (1.0 - pp * pp).sqrt();
            let sdel = f64::from((j as f32 * RAD).sin());
            let gg = v[l] * sh * dtih / v[1] / sdel / cs0;
            spreading[j] = gg.sqrt() / 6371.0 * 100000.0;
        }
        (spreading, slowness)
    }
}

/// The distance in radians of a ray with parameter `p` from a source in
/// layer `l` at radius `rh`: down to its turning point, twice, then back up
/// past the source to the surface.
fn ray_distance(p: f64, l: usize, rh: f32, v: &[f64], r: &[f32]) -> f64 {
    let mut del = 0.0;
    let (mut r1, mut r2, mut i) = (rh, r[l], l);
    while p * v[i] / f64::from(r1) < 1.0 {
        let del1 = (p * v[i] / f64::from(r1)).asin();
        let pa = p * v[i] / f64::from(r2);
        if pa >= 1.0 {
            del = del + f64::from(PI / 2.0) - del1;
            break;
        }
        del = del + pa.asin() - del1;
        r1 = r[i];
        i += 1;
        r2 = r[i];
    }
    del *= 2.0;
    let mut r2 = rh;
    for i in (1..l).rev() {
        let r1 = r[i];
        del = del + (p * v[i] / f64::from(r2)).asin() - (p * v[i] / f64::from(r1)).asin();
        r2 = r1;
    }
    del + (p * v[1] / 6371.0).asin() - (p * v[1] / f64::from(r2)).asin()
}

/// The tables for one source depth.
pub struct Earth<'t> {
    tables: &'t Tables,
    depth: f32,
    /// Indexed by whole degree.
    spreading: Vec<f64>,
    slowness: Vec<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Geometry {
    /// Station azimuth from the source, in degrees.
    pub azimuth: f32,
    /// Source azimuth from the station, in degrees.
    pub back_azimuth: f32,
    /// Geocentric distance by Utsu's formula, in degrees.
    pub distance: f32,
    /// Ray parameter in s/km, at the spherical distance.
    pub slowness: f32,
    /// Geometrical spreading, at the spherical distance.
    pub spreading: f32,
}

impl Earth<'_> {
    pub fn depth(&self) -> f32 {
        self.depth
    }

    /// `parm`: the geometry from a source to a station, each given as
    /// (latitude, longitude) in degrees.
    pub fn geometry(&self, source: (f32, f32), station: (f32, f32)) -> Geometry {
        let t0 = source.0 * RAD;
        let t1 = station.0 * RAD;
        let dphi = (station.1 - source.1) * RAD;
        let (c0, c1, s0, s1) = (t0.cos(), t1.cos(), t0.sin(), t1.sin());
        let (c2, s2) = (dphi.cos(), dphi.sin());
        let sphere = (c0 * c1 * c2 + s0 * s1).acos() / RAD;
        // Utsu's correction from geographic to geocentric latitude.
        const UTSU: f32 = 11.55 / 60.0 * RAD;
        let t00 = t0 - UTSU * (2.0 * t0).sin();
        let t11 = t1 - UTSU * (2.0 * t1).sin();
        let d11 = t00.cos() * t11.cos() * c2 + t00.sin() * t11.sin();
        let distance = d11.acos() / RAD;
        let azimuth = (c1 * s2).atan2(s1 * c0 - s0 * c1 * c2) / RAD;
        let back_azimuth = (-c0 * s2).atan2(s0 * c1 - s1 * c0 * c2) / RAD;
        let whole = sphere as usize;
        let fraction = f64::from(sphere - whole as f32);
        let at = |table: &[f64]| {
            let (a, b) = (table[whole], table[(whole + 1).min(SPAN - 1)]);
            (a + fraction * (b - a)) as f32
        };
        Geometry {
            azimuth,
            back_azimuth,
            distance,
            slowness: at(&self.slowness),
            spreading: at(&self.spreading),
        }
    }

    /// `jbtravel`: seconds from origin to the phase onset at `distance`
    /// degrees. The tables cover 10 to 103 degrees; outside, it is 0.
    pub fn travel(&self, phase: Phase, distance: f32) -> f32 {
        if !(10.0..=103.0).contains(&distance) {
            return 0.0;
        }
        let table = match phase {
            Phase::P => &self.tables.p,
            Phase::S => &self.tables.s,
        };
        let c = DEPTHS.iter().position(|&h| self.depth < h).unwrap_or(13) - 1;
        let ratio = (self.depth - DEPTHS[c]) / (DEPTHS[c + 1] - DEPTHS[c]);
        let whole = distance as usize;
        let t1 = table.at(whole, c, ratio);
        let t2 = table.at(whole + 1, c, ratio);
        t1 + (distance - whole as f32) * (t2 - t1)
    }
}
