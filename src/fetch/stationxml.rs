//! The channels of an FDSN StationXML document at `level=response`, with the
//! poles and zeros `prep` divides out.

use roxmltree::{Document, Node};

use crate::prep::fft::Complex;
use crate::prep::{Quantity, Response};

pub struct Channel {
    pub network: String,
    pub station: String,
    pub location: String,
    pub code: String,
    pub latitude: f32,
    pub longitude: f32,
    /// The response, or why it cannot be used.
    pub response: Result<Response, String>,
}

impl Channel {
    /// `NET.STA.LOC.CHA`.
    pub fn id(&self) -> String {
        format!(
            "{}.{}.{}.{}",
            self.network, self.station, self.location, self.code
        )
    }
}

/// Every channel epoch in the document, in document order.
pub fn parse(xml: &str) -> Result<Vec<Channel>, String> {
    let doc = Document::parse(xml).map_err(|e| e.to_string())?;
    let mut channels = Vec::new();
    for network in children(doc.root_element(), "Network") {
        for station in children(network, "Station") {
            for channel in children(station, "Channel") {
                channels.push(Channel {
                    network: attribute(network, "code")?,
                    station: attribute(station, "code")?,
                    location: attribute(channel, "locationCode")?,
                    code: attribute(channel, "code")?,
                    latitude: number(channel, "Latitude")?,
                    longitude: number(channel, "Longitude")?,
                    response: child(channel, "Response").and_then(response),
                });
            }
        }
    }
    Ok(channels)
}

/// The first poles-and-zeros stage in rad/s, normalised as the overall
/// sensitivity expects: counts = sensitivity * normalisation * prod(s - z) /
/// prod(s - p).
fn response(node: Node) -> Result<Response, String> {
    let sensitivity = child(node, "InstrumentSensitivity")?;
    let stage = children(node, "Stage")
        .find_map(|s| child(s, "PolesZeros").ok())
        .ok_or("no poles-and-zeros stage")?;
    // A stage in Hz has roots and a variable 2 pi times smaller; scaling the
    // roots by 2 pi scales the transfer function by 2 pi ^ (zeros - poles).
    let hertz = match text(stage, "PzTransferFunctionType")?.as_str() {
        "LAPLACE (RADIANS/SECOND)" => 1.0,
        "LAPLACE (HERTZ)" => 2.0 * std::f64::consts::PI,
        kind => return Err(format!("transfer function {kind} is not Laplace")),
    };
    let roots = |name| {
        children(stage, name)
            .map(|n| {
                let re = f64::from(number(n, "Real")?) * hertz;
                let im = f64::from(number(n, "Imaginary")?) * hertz;
                Ok(Complex::new(re as f32, im as f32))
            })
            .collect::<Result<Vec<_>, String>>()
    };
    let poles = roots("Pole")?;
    let zeros = roots("Zero")?;
    let excess = poles.len() as i32 - zeros.len() as i32;
    let normalisation = f64::from(number(stage, "NormalizationFactor")?) * hertz.powi(excess);
    Ok(Response {
        units: units(&text(child(stage, "InputUnits")?, "Name")?)?,
        normalisation: normalisation as f32,
        sensitivity: number(sensitivity, "Value")?,
        poles,
        zeros,
    })
}

fn units(name: &str) -> Result<Quantity, String> {
    match name.to_ascii_lowercase().as_str() {
        "m" => Ok(Quantity::Displacement),
        "m/s" => Ok(Quantity::Velocity),
        "m/s**2" | "m/s/s" | "m/s2" => Ok(Quantity::Acceleration),
        _ => Err(format!("input units {name:?} are not m, m/s or m/s**2")),
    }
}

fn children<'a, 'i>(node: Node<'a, 'i>, name: &'static str) -> impl Iterator<Item = Node<'a, 'i>> {
    node.children().filter(move |n| n.has_tag_name(name))
}

fn child<'a, 'i>(node: Node<'a, 'i>, name: &str) -> Result<Node<'a, 'i>, String> {
    node.children()
        .find(|n| n.has_tag_name(name))
        .ok_or_else(|| format!("<{}> has no <{name}>", node.tag_name().name()))
}

fn text(node: Node, name: &str) -> Result<String, String> {
    Ok(child(node, name)?.text().unwrap_or("").trim().to_owned())
}

fn number(node: Node, name: &str) -> Result<f32, String> {
    let value = text(node, name)?;
    value
        .parse()
        .map_err(|_| format!("<{name}> {value:?} is not a number"))
}

fn attribute(node: Node, name: &str) -> Result<String, String> {
    node.attribute(name)
        .map(str::to_owned)
        .ok_or_else(|| format!("<{}> has no {name}", node.tag_name().name()))
}
