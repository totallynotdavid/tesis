//! What `inversion` writes: the `fort.16` report, the `fort.4` plot data
//! and its one line on standard output.

use super::fortran::{f, i};
use super::inversion::{Basis, TimeFunction};
use super::unformatted::{padded, reals, record};

#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    /// The `i_inversion` title.
    pub title: String,
    /// The `fort.3` title, 40 characters.
    pub structure: String,
    pub dt: f32,
    /// Depth of the reference point, in km.
    pub h0: f32,
    /// Dip of the line of depths, in degrees.
    pub dip: f32,
    pub basis: Basis,
    /// Samples added to every onset.
    pub shift: i32,
    pub time_function: TimeFunction,
    pub max_velocity: f32,
    /// The latest onset searched for each subevent, in seconds.
    pub spans: Vec<f32>,
    /// The source time function, one value per fitted sample.
    pub source: Vec<f32>,
    /// The line of depths: spacing in km, count, and the point at `h0`.
    pub dk: f32,
    pub nk: usize,
    pub k0: i32,
    /// The grid points at each depth, as (distance, azimuth).
    pub points: Vec<(f32, f32)>,
    pub stations: Vec<Fit>,
    /// The subevents as each round of the search left them; the last round
    /// is the result.
    pub rounds: Vec<Vec<Subevent>>,
    /// The residual variance after each subevent of the first round,
    /// relative to the data.
    pub convergence: Vec<f32>,
    /// The final residual variance, relative to the data.
    pub variance: f32,
    /// The sum of the subevents.
    pub total: Mechanism,
    /// The change in `total` when each weight moves by one standard error.
    pub deviation: Mechanism,
    /// The summed tensor, nine reals in Fortran column order.
    pub tensor: [f32; 9],
    /// The tension axis of the summed tensor.
    pub axis: [f32; 3],
    pub magnitude: f32,
}

/// One station the inversion fitted.
#[derive(Debug, Clone, PartialEq)]
pub struct Fit {
    /// 12 characters.
    pub name: String,
    /// `im`, `ib` and `ic` from `fort.1`.
    pub instrument: i32,
    pub phase: i32,
    pub component: i32,
    pub azimuth: f32,
    pub back_azimuth: f32,
    pub distance: f32,
    pub slowness: f32,
    pub spreading: f32,
    pub weight: f32,
    /// The onset with the shift added.
    pub onset: i32,
    /// The source-layer velocity of the station's phase.
    pub velocity: f32,
    /// The record and the final residual, in micrometres at a common gain.
    pub observed: Vec<f32>,
    pub residual: Vec<f32>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Subevent {
    /// Onset in seconds.
    pub time: f32,
    /// The grid point: distance and azimuth from the hypocentre, and the
    /// offset along the line of depths, in km.
    pub distance: f32,
    pub azimuth: f32,
    pub offset: f32,
    pub mechanism: Mechanism,
    /// The residual variance once this subevent is removed.
    pub variance: f32,
    /// The tensor, nine reals in Fortran column order.
    pub tensor: [f32; 9],
    pub axis: [f32; 3],
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mechanism {
    /// In units of 1e18 N m.
    pub moment: f32,
    /// The non-double-couple part; the isotropic weight under
    /// [`Basis::MomentTensor`].
    pub second: f32,
    pub strike: f32,
    pub dip: f32,
    pub slip: f32,
}

/// The `fort.16` report.
pub fn report(s: &Solution) -> String {
    let mut out = String::new();
    let mut line = |text: String| {
        out.push_str(&text);
        out.push('\n');
    };
    line(format!(" {:<40.40}", s.title));
    line(format!(" {:<40.40}", s.structure));
    line("   NT    DT     H0    Dip    Ne Type Shift".into());
    line(format!(
        "{}{}{}{}{}{}{}",
        i(s.source.len() as i32, 6),
        f(s.dt, 6, 3),
        f(s.h0, 7, 1),
        f(s.dip, 7, 1),
        i(s.spans.len() as i32, 5),
        i(s.basis.code(), 5),
        i(s.shift, 5)
    ));
    let t = s.time_function;
    line(format!(
        " Time Function : (Type,T1,T2) = {}{}{}",
        i(t.shape.code(), 1),
        f(t.t1, 6, 1),
        f(t.t2, 6, 1)
    ));
    for text in extent(s.max_velocity, &s.spans) {
        line(text);
    }
    let second = if s.basis == Basis::MomentTensor {
        "Isotrpc"
    } else {
        "non.DC "
    };
    line(format!(
        "  # Time Location(r,phi,Y-axis)  Moment {second} Strike  Dip   Slip  Error"
    ));
    for (round, subevents) in s.rounds.iter().enumerate() {
        line(format!("Step : {}", i(round as i32 + 1, 2)));
        for (n, e) in subevents.iter().enumerate() {
            let m = e.mechanism;
            line(format!(
                "{}{}{}{}{}{}{}{}{}{}{}",
                i(n as i32 + 1, 4),
                f(e.time, 6, 1),
                f(e.distance, 7, 1),
                f(e.azimuth, 7, 1),
                f(e.offset, 7, 1),
                f(m.moment, 9, 3),
                f(m.second, 9, 3),
                f(m.strike, 7, 1),
                f(m.dip, 7, 1),
                f(m.slip, 7, 1),
                f(e.variance, 8, 4)
            ));
        }
    }
    line("< Convergence > ".into());
    for (n, &e) in s.convergence.iter().enumerate() {
        let bar = (e * 50.0 + 1.0).round() as i32;
        let bar = if bar > 0 {
            format!("  {}", "=".repeat(bar as usize))
        } else {
            String::new()
        };
        line(format!("     {}{}{bar}", i(n as i32 + 1, 4), f(e, 8, 4)));
    }
    let row = |label: &str, m: Mechanism| {
        format!(
            "          {label}{}{}{}{}{}",
            f(m.moment, 8, 2),
            f(m.second, 8, 2),
            f(m.strike, 8, 1),
            f(m.dip, 8, 1),
            f(m.slip, 8, 1)
        )
    };
    line(row("Total                ", s.total));
    line(row("standard-deviation   ", s.deviation));
    line(format!(" Mw = {}", f(s.magnitude, 7, 2)));
    line(" Stn          Inst Mode Comp   Az     B.Az   Del    p      G      fc    shift".into());
    for st in &s.stations {
        line(format!(
            "{:<12.12}{}{}{}{}{}{}{}{}{}{}",
            st.name,
            i(st.instrument, 5),
            i(st.phase, 5),
            i(st.component, 5),
            f(st.azimuth, 7, 1),
            f(st.back_azimuth, 7, 1),
            f(st.distance, 7, 1),
            f(st.slowness, 8, 3),
            f(st.spreading, 7, 2),
            f(st.weight, 7, 2),
            i(st.onset, 6)
        ));
    }
    out
}

/// Format 203, which reverts to its start after twelve spans: a velocity
/// line, then a span line, again for each further value, stopping at the
/// first descriptor with no value left.
fn extent(max_velocity: f32, spans: &[f32]) -> Vec<String> {
    let mut values = std::iter::once(max_velocity).chain(spans.iter().copied());
    let mut lines = Vec::new();
    let mut line = String::new();
    'format: loop {
        line.push_str(" Grid-search extent: V max = ");
        let Some(v) = values.next() else { break };
        line.push_str(&f(v, 6, 1));
        lines.push(std::mem::take(&mut line));
        line.push_str(" Time span at each iteration = ");
        for _ in 0..12 {
            let Some(v) = values.next() else {
                break 'format;
            };
            line.push_str(&f(v, 6, 1));
        }
        let mut rest = values.clone().peekable();
        if rest.peek().is_none() {
            break;
        }
        lines.push(std::mem::take(&mut line));
    }
    lines.push(line);
    lines
}

/// The `fort.4` plot data.
pub fn plot(s: &Solution) -> Vec<u8> {
    let mut out = Vec::new();
    let ints =
        |values: &[i32]| -> Vec<u8> { values.iter().flat_map(|v| v.to_le_bytes()).collect() };
    let t = s.time_function;
    record(&mut out, &padded(&s.title, 40));
    let mut header = ints(&[s.stations.len() as i32, s.source.len() as i32]);
    header.extend(reals(&[s.dt]));
    header.extend(ints(&[s.spans.len() as i32]));
    header.extend(reals(&[t.t1, t.t2]));
    header.extend(ints(&[t.shape.code()]));
    record(&mut out, &header);
    record(&mut out, &reals(&s.source));
    let mut line = ints(&[s.points.len() as i32]);
    line.extend(reals(&[s.dk]));
    line.extend(ints(&[s.nk as i32, s.k0]));
    record(&mut out, &line);
    let points: Vec<f32> = s.points.iter().flat_map(|&(r, phi)| [r, phi]).collect();
    record(&mut out, &reals(&points));
    for st in &s.stations {
        let mut name = padded(&st.name, 12);
        name.extend(ints(&[st.instrument, st.phase, st.component]));
        record(&mut out, &name);
        record(&mut out, &reals(&[st.azimuth, st.slowness * st.velocity]));
        record(&mut out, &reals(&st.observed));
    }
    let last = s.rounds.last().map(Vec::as_slice).unwrap_or(&[]);
    for (n, e) in last.iter().enumerate() {
        let m = e.mechanism;
        let mut row = ints(&[n as i32 + 1]);
        row.extend(reals(&[
            e.time, e.distance, e.azimuth, e.offset, m.moment, m.second, m.strike, m.dip, m.slip,
        ]));
        record(&mut out, &row);
        record(&mut out, &reals(&e.tensor));
        record(&mut out, &reals(&e.axis));
    }
    record(&mut out, &reals(&[s.variance]));
    for st in &s.stations {
        record(&mut out, &reals(&st.residual));
    }
    let m = s.total;
    record(
        &mut out,
        &reals(&[m.moment, m.second, m.strike, m.dip, m.slip]),
    );
    record(&mut out, &reals(&s.tensor));
    record(&mut out, &reals(&s.axis));
    out
}

/// The line `inversion` prints on standard output.
pub fn summary(s: &Solution) -> String {
    format!(
        " Mw, Var. = {}{}\n",
        f(s.magnitude, 7, 2),
        f(s.variance, 8, 4)
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_extent_reverts_to_a_velocity_line_after_twelve_spans() {
        assert_eq!(
            extent(4.0, &[100.0; 2]),
            [
                " Grid-search extent: V max =    4.0",
                " Time span at each iteration =  100.0 100.0"
            ]
        );
        let lines = extent(4.0, &[1.0; 13]);
        assert_eq!(lines.len(), 4);
        assert_eq!(lines[2], " Grid-search extent: V max =    1.0");
        assert_eq!(lines[3], " Time span at each iteration = ");
    }
}
