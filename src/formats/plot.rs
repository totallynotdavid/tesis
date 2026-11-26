//! `fort.4`: the plot data `inversion` writes and `graphics` reads, a
//! gfortran unformatted file.

use super::inversion::{Shape, TimeFunction};
use super::solution::{Mechanism, Solution};
use super::unformatted::{Error, Reader, padded, reals, record, word};

const TITLE: usize = 40;
const NAME: usize = 12;

#[derive(Debug, Clone, PartialEq)]
pub struct Plot {
    /// The `i_inversion` title, 40 characters.
    pub title: String,
    pub dt: f32,
    pub time_function: TimeFunction,
    /// The source time function, one value per fitted sample.
    pub source: Vec<f32>,
    /// The line of depths: spacing in km, count, and the point at `h0`.
    pub dk: f32,
    pub nk: usize,
    pub k0: i32,
    /// The grid points at each depth, as (distance, azimuth).
    pub points: Vec<(f32, f32)>,
    pub stations: Vec<Trace>,
    pub subevents: Vec<Subevent>,
    /// The final residual variance, relative to the data.
    pub variance: f32,
    pub total: Mechanism,
    /// The summed tensor, nine reals in Fortran column order.
    pub tensor: [f32; 9],
    /// The tension axis of the summed tensor.
    pub axis: [f32; 3],
}

/// One fitted station.
#[derive(Debug, Clone, PartialEq)]
pub struct Trace {
    /// 12 characters.
    pub name: String,
    /// `im`, `ib` and `ic` from `fort.1`.
    pub instrument: i32,
    pub phase: i32,
    pub component: i32,
    pub azimuth: f32,
    /// The sine of the take-off angle at the source.
    pub ray: f32,
    /// The record and the final residual, `source.len()` samples each.
    pub observed: Vec<f32>,
    pub residual: Vec<f32>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Subevent {
    /// Onset in seconds.
    pub time: f32,
    /// Distance and azimuth from the hypocentre, and the offset along the
    /// line of depths, in km.
    pub distance: f32,
    pub azimuth: f32,
    pub offset: f32,
    pub mechanism: Mechanism,
    /// The tensor, nine reals in Fortran column order.
    pub tensor: [f32; 9],
    pub axis: [f32; 3],
}

impl From<&Solution> for Plot {
    /// The last round's subevents, as `inversion` writes them.
    fn from(s: &Solution) -> Self {
        let last = s.rounds.last().map(Vec::as_slice).unwrap_or(&[]);
        Self {
            title: s.title.clone(),
            dt: s.dt,
            time_function: s.time_function,
            source: s.source.clone(),
            dk: s.dk,
            nk: s.nk,
            k0: s.k0,
            points: s.points.clone(),
            stations: s
                .stations
                .iter()
                .map(|st| Trace {
                    name: st.name.clone(),
                    instrument: st.instrument,
                    phase: st.phase,
                    component: st.component,
                    azimuth: st.azimuth,
                    ray: st.slowness * st.velocity,
                    observed: st.observed.clone(),
                    residual: st.residual.clone(),
                })
                .collect(),
            subevents: last
                .iter()
                .map(|e| Subevent {
                    time: e.time,
                    distance: e.distance,
                    azimuth: e.azimuth,
                    offset: e.offset,
                    mechanism: e.mechanism,
                    tensor: e.tensor,
                    axis: e.axis,
                })
                .collect(),
            variance: s.variance,
            total: s.total,
            tensor: s.tensor,
            axis: s.axis,
        }
    }
}

pub fn write(p: &Plot) -> Vec<u8> {
    let mut out = Vec::new();
    let t = p.time_function;
    record(&mut out, &padded(&p.title, TITLE));
    let mut header = ints(&[p.stations.len() as i32, p.source.len() as i32]);
    header.extend(reals(&[p.dt]));
    header.extend(ints(&[p.subevents.len() as i32]));
    header.extend(reals(&[t.t1, t.t2]));
    header.extend(ints(&[t.shape.code()]));
    record(&mut out, &header);
    record(&mut out, &reals(&p.source));
    let mut line = ints(&[p.points.len() as i32]);
    line.extend(reals(&[p.dk]));
    line.extend(ints(&[p.nk as i32, p.k0]));
    record(&mut out, &line);
    let points: Vec<f32> = p.points.iter().flat_map(|&(r, phi)| [r, phi]).collect();
    record(&mut out, &reals(&points));
    for st in &p.stations {
        let mut name = padded(&st.name, NAME);
        name.extend(ints(&[st.instrument, st.phase, st.component]));
        record(&mut out, &name);
        record(&mut out, &reals(&[st.azimuth, st.ray]));
        record(&mut out, &reals(&st.observed));
    }
    for (n, e) in p.subevents.iter().enumerate() {
        let mut row = ints(&[n as i32 + 1]);
        row.extend(reals(&[e.time, e.distance, e.azimuth, e.offset]));
        row.extend(reals(&mechanism(e.mechanism)));
        record(&mut out, &row);
        record(&mut out, &reals(&e.tensor));
        record(&mut out, &reals(&e.axis));
    }
    record(&mut out, &reals(&[p.variance]));
    for st in &p.stations {
        record(&mut out, &reals(&st.residual));
    }
    record(&mut out, &reals(&mechanism(p.total)));
    record(&mut out, &reals(&p.tensor));
    record(&mut out, &reals(&p.axis));
    out
}

pub fn read(bytes: &[u8]) -> Result<Plot, Error> {
    let mut records = Reader::new(bytes);
    let title = records.next_text(TITLE)?.trim_end().to_owned();
    let header = records.next(28)?;
    let (int, real) = (int_at(&header), real_at(&header));
    let js = records.count(int(0))?;
    let nt = records.count(int(1))?;
    let ite = records.count(int(3))?;
    let shape = Shape::from_code(int(6))
        .ok_or_else(|| records.error(format!("time function {} is not 0 to 3", int(6))))?;
    let time_function = TimeFunction {
        shape,
        t1: real(4),
        t2: real(5),
    };
    let source = records.next_reals(nt)?;
    let line = records.next(16)?;
    let (int, real) = (int_at(&line), real_at(&line));
    let nlen = records.count(int(0))?;
    let (dk, nk, k0) = (real(1), records.count(int(2))?, int(3));
    let points = records.next_reals(2 * nlen)?;
    let points = points.chunks(2).map(|p| (p[0], p[1])).collect();
    let mut stations = Vec::with_capacity(js);
    for _ in 0..js {
        let head = records.next(NAME + 12)?;
        let int = int_at(&head[NAME..]);
        let geometry = records.next_reals(2)?;
        stations.push(Trace {
            name: String::from_utf8_lossy(&head[..NAME]).into_owned(),
            instrument: int(0),
            phase: int(1),
            component: int(2),
            azimuth: geometry[0],
            ray: geometry[1],
            observed: records.next_reals(nt)?,
            residual: Vec::new(),
        });
    }
    let mut subevents = Vec::with_capacity(ite);
    for _ in 0..ite {
        let row = records.next(40)?;
        let real = real_at(&row);
        subevents.push(Subevent {
            time: real(1),
            distance: real(2),
            azimuth: real(3),
            offset: real(4),
            mechanism: Mechanism {
                moment: real(5),
                second: real(6),
                strike: real(7),
                dip: real(8),
                slip: real(9),
            },
            tensor: array(records.next_reals(9)?),
            axis: array(records.next_reals(3)?),
        });
    }
    let variance = records.next_reals(1)?[0];
    for st in &mut stations {
        st.residual = records.next_reals(nt)?;
    }
    let total = records.next_reals(5)?;
    Ok(Plot {
        title,
        dt: real_at(&header)(2),
        time_function,
        source,
        dk,
        nk,
        k0,
        points,
        stations,
        subevents,
        variance,
        total: Mechanism {
            moment: total[0],
            second: total[1],
            strike: total[2],
            dip: total[3],
            slip: total[4],
        },
        tensor: array(records.next_reals(9)?),
        axis: array(records.next_reals(3)?),
    })
}

fn ints(values: &[i32]) -> Vec<u8> {
    values.iter().flat_map(|v| v.to_le_bytes()).collect()
}

fn mechanism(m: Mechanism) -> [f32; 5] {
    [m.moment, m.second, m.strike, m.dip, m.slip]
}

fn int_at(bytes: &[u8]) -> impl Fn(usize) -> i32 + '_ {
    |i| i32::from_le_bytes(word(&bytes[4 * i..]))
}

fn real_at(bytes: &[u8]) -> impl Fn(usize) -> f32 + '_ {
    |i| f32::from_le_bytes(word(&bytes[4 * i..]))
}

/// A record of exactly `N` reals, which `next_reals(N)` guarantees.
fn array<const N: usize>(values: Vec<f32>) -> [f32; N] {
    values.try_into().unwrap()
}
