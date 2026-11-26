//! `i_inversion`: the parameters of `inversion`, read from its standard
//! input. The example file holds a second parameter set after the first; the
//! program reads only the first.

use std::fmt::Write as _;

use serde::Deserialize;

use super::reals;
use super::records::{Error, Records};

#[derive(Debug, Clone, PartialEq)]
pub struct Parameters {
    /// The `fort.16` title, at most 40 characters.
    pub title: String,
    /// Length of the fitted waveforms in seconds.
    pub length: f32,
    /// Subevents start no earlier than their distance from the hypocentre
    /// over this velocity, in km/s.
    pub max_velocity: f32,
    pub time_function: TimeFunction,
    /// For each subevent, the latest onset searched, in seconds.
    pub spans: Vec<f32>,
    /// One weight per station `green` kept, in `fort.1` order. A station
    /// past the list, or with weight 0, is left out.
    pub weights: Vec<f32>,
    /// Samples added to every station's onset.
    pub shift: i32,
    pub basis: Basis,
    pub grid: Grid,
}

/// The source time function of each subevent.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TimeFunction {
    pub shape: Shape,
    #[serde(default)]
    pub t1: f32,
    #[serde(default)]
    pub t2: f32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Shape {
    /// One sample of area 1.
    Impulse,
    /// Rise `t1`, top `t2 - t1`, fall `t1`; a ramp of rise `t1` if `t2 < t1`.
    Trapezoid,
    /// The trapezoid with cosine flanks.
    Cosine,
    /// A triangle peaking at `t1` and ending at `t2`.
    Triangle,
}

/// The moment-tensor elements a subevent combines.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Basis {
    /// The five deviatoric elements.
    Deviatoric,
    /// The five deviatoric elements, constrained to a double couple.
    DoubleCouple,
    /// Mxy and Mxx-Myy: vertical strike-slip.
    StrikeSlip,
    /// All six elements, with the isotropic part.
    MomentTensor,
    /// Four elements, constrained to a double couple with a vertical nodal
    /// plane.
    VerticalDoubleCouple,
}

/// The source points searched at each depth.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields, rename_all = "lowercase")]
pub enum Grid {
    /// Points given as (distance in km, azimuth in degrees) from the
    /// hypocentre.
    Points(Vec<(f32, f32)>),
    /// `count` points along the azimuth `azimuth`, `spacing` km apart, with
    /// point `reference`, counted from 1, at the hypocentre.
    Line {
        azimuth: f32,
        spacing: f32,
        count: usize,
        #[serde(default = "first")]
        reference: i32,
    },
}

fn first() -> i32 {
    1
}

impl Shape {
    const CODES: [Self; 4] = [Self::Impulse, Self::Trapezoid, Self::Cosine, Self::Triangle];

    pub fn code(self) -> i32 {
        Self::CODES.iter().position(|&s| s == self).unwrap() as i32
    }

    fn from_code(code: i32) -> Option<Self> {
        Self::CODES.get(usize::try_from(code).ok()?).copied()
    }
}

impl Basis {
    const CODES: [Self; 5] = [
        Self::Deviatoric,
        Self::DoubleCouple,
        Self::StrikeSlip,
        Self::MomentTensor,
        Self::VerticalDoubleCouple,
    ];

    /// `isr` in the reference.
    pub fn code(self) -> i32 {
        Self::CODES.iter().position(|&b| b == self).unwrap() as i32
    }

    fn from_code(code: i32) -> Option<Self> {
        Self::CODES.get(usize::try_from(code).ok()?).copied()
    }

    /// How many of the base elements, in `fort.3` order, it uses.
    pub fn elements(self) -> usize {
        match self {
            Self::Deviatoric | Self::DoubleCouple => 5,
            Self::StrikeSlip => 2,
            Self::MomentTensor => 6,
            Self::VerticalDoubleCouple => 4,
        }
    }

    pub fn is_double_couple(self) -> bool {
        matches!(self, Self::DoubleCouple | Self::VerticalDoubleCouple)
    }
}

impl Grid {
    /// The points as (distance, azimuth).
    pub fn points(&self) -> Vec<(f32, f32)> {
        match self {
            Self::Points(points) => points.clone(),
            &Self::Line {
                azimuth,
                spacing,
                count,
                reference,
            } => (1..=count as i32)
                .map(|l| (spacing * (l - reference) as f32, azimuth))
                .collect(),
        }
    }
}

pub fn parse(text: &str) -> Result<Parameters, Error> {
    let mut records = Records::new(text);
    let line = records.raw().unwrap_or_default();
    let title: String = line
        .chars()
        .take(40)
        .collect::<String>()
        .trim_end()
        .to_owned();
    let mut read = records.read();
    let (length, max_velocity) = (read.f32()?, read.f32()?);
    drop(read);
    let mut read = records.read();
    let (shape, t1, t2) = (read.i32()?, read.f32()?, read.f32()?);
    drop(read);
    let shape = Shape::from_code(shape)
        .ok_or_else(|| records.error(format!("time function {shape} is not 0 to 3")))?;
    let time_function = TimeFunction { shape, t1, t2 };
    let subevents = records.read().count()?;
    let spans = records.read().f32s(subevents)?;
    let mut read = records.read();
    let stations = read.count()?;
    let weights = read.f32s(stations)?;
    drop(read);
    let mut read = records.read();
    let (shift, isr, ifp) = (read.i32()?, read.i32()?, read.i32()?);
    drop(read);
    let basis =
        Basis::from_code(isr).ok_or_else(|| records.error(format!("basis {isr} is not 0 to 4")))?;
    let grid = if ifp == 0 {
        let count = records.read().count()?;
        let distances = records.read().f32s(count)?;
        let azimuths = records.read().f32s(count)?;
        Grid::Points(distances.into_iter().zip(azimuths).collect())
    } else {
        let mut read = records.read();
        Grid::Line {
            azimuth: read.f32()?,
            spacing: read.f32()?,
            count: read.count()?,
            reference: read.i32()?,
        }
    };
    Ok(Parameters {
        title,
        length,
        max_velocity,
        time_function,
        spans,
        weights,
        shift,
        basis,
        grid,
    })
}

pub fn write(parameters: &Parameters) -> String {
    let p = parameters;
    let f = p.time_function;
    let mut out = format!("{}\n", p.title);
    let _ = writeln!(out, "{}", reals(&[p.length, p.max_velocity]));
    let _ = writeln!(out, "{}{}", f.shape.code(), reals(&[f.t1, f.t2]));
    let _ = writeln!(out, "{}", p.spans.len());
    let _ = writeln!(out, "{}", reals(&p.spans));
    let _ = writeln!(out, "{}{}", p.weights.len(), reals(&p.weights));
    match &p.grid {
        Grid::Points(points) => {
            let _ = writeln!(out, "{} {} 0", p.shift, p.basis.code());
            let _ = writeln!(out, "{}", points.len());
            let (distances, azimuths): (Vec<f32>, Vec<f32>) = points.iter().copied().unzip();
            let _ = writeln!(out, "{}", reals(&distances));
            let _ = writeln!(out, "{}", reals(&azimuths));
        }
        &Grid::Line {
            azimuth,
            spacing,
            count,
            reference,
        } => {
            let _ = writeln!(out, "{} {} 1", p.shift, p.basis.code());
            let _ = writeln!(out, "{} {count} {reference}", reals(&[azimuth, spacing]));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_the_first_parameter_set_of_the_example() {
        let text = include_str!("../../reference/examples/input/i_inversion");
        let parameters = parse(text).unwrap();
        assert_eq!(parameters.spans, [100.0; 3]);
        assert_eq!(parameters.weights.len(), 24);
        assert_eq!(parameters.weights[23], 0.1);
        assert_eq!(parameters.shift, -30);
        assert_eq!(parameters.basis, Basis::DoubleCouple);
        assert_eq!(parameters.grid.points()[1], (30.0, 120.0));
        assert_eq!(parse(&write(&parameters)).unwrap(), parameters);
    }
}
