//! `fort.2`: the layered velocity structure at the source, under the
//! receivers, and at the PP bounce point. The programs read the first model in
//! the file; the example `fort.2` holds several more after it.

use std::fmt::Write as _;

use super::reals;
use super::records::{Error, Read, Records};

#[derive(Debug, Clone, PartialEq)]
pub struct Structure {
    pub title: String,
    /// Attenuation t* for P and S, in seconds.
    pub tq_p: f32,
    pub tq_s: f32,
    pub source: Vec<Layer>,
    pub receiver: Vec<Layer>,
    pub bounce: Vec<Layer>,
}

/// One layer, top down. The last layer is the half-space; its thickness is
/// ignored.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Layer {
    /// P and S velocity in km/s.
    pub vp: f32,
    pub vs: f32,
    /// Density in g/cm3.
    pub density: f32,
    /// Thickness in km.
    pub thickness: f32,
}

pub fn parse(text: &str) -> Result<Structure, Error> {
    let mut records = Records::new(text);
    let title = records.raw().unwrap_or_default().trim_end().to_owned();
    let mut read = records.read();
    let (tq_p, tq_s) = (read.f32()?, read.f32()?);
    let source = layers(&mut read)?;
    drop(read);
    let receiver = layers(&mut records.read())?;
    let bounce = layers(&mut records.read())?;
    Ok(Structure {
        title,
        tq_p,
        tq_s,
        source,
        receiver,
        bounce,
    })
}

fn layers(read: &mut Read) -> Result<Vec<Layer>, Error> {
    let n = read.count()?;
    (0..n)
        .map(|_| {
            Ok(Layer {
                vp: read.f32()?,
                vs: read.f32()?,
                density: read.f32()?,
                thickness: read.f32()?,
            })
        })
        .collect()
}

pub fn write(structure: &Structure) -> String {
    let mut out = format!("{}\n", structure.title);
    let head = reals(&[structure.tq_p, structure.tq_s]);
    write_layers(&mut out, &head, &structure.source);
    write_layers(&mut out, "", &structure.receiver);
    write_layers(&mut out, "", &structure.bounce);
    out
}

/// One layer per record, the count before the first. Fortran reads the
/// whole group as one statement, so the line breaks carry no meaning.
fn write_layers(out: &mut String, head: &str, layers: &[Layer]) {
    let _ = write!(out, "{head} {:>3}", layers.len());
    for layer in layers {
        let _ = writeln!(
            out,
            "{}",
            reals(&[layer.vp, layer.vs, layer.density, layer.thickness])
        );
    }
    if layers.is_empty() {
        out.push('\n');
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TWO_MODELS: &str = "J-B model
1.0 4.0 2 5.57 3.36 2.65 15.0
          8.10 4.68 3.30 0.0
        1 6.50 3.74 2.87 18.0
        1 6.50 3.74 2.87 0.0

second model
1.0 4.0 1 5.57 3.36 2.65 15.0
";

    #[test]
    fn only_the_first_model_is_read() {
        let structure = parse(TWO_MODELS).unwrap();
        assert_eq!(structure.source.len(), 2);
        assert_eq!(structure.source[1].vp, 8.1);
        assert_eq!(structure.bounce[0].thickness, 0.0);
    }

    #[test]
    fn writing_then_reading_gives_the_same_values() {
        let structure = parse(TWO_MODELS).unwrap();
        assert_eq!(parse(&write(&structure)).unwrap(), structure);
    }
}
