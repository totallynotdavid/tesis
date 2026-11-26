//! A minimal SVG writer. Coordinates are pixels, y down.

use std::fmt::Write as _;

pub struct Svg {
    body: String,
}

impl Svg {
    pub fn new() -> Self {
        Self {
            body: String::new(),
        }
    }

    /// The document, `width` by `height` pixels on a white page.
    pub fn finish(self, width: f64, height: f64) -> String {
        format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{h}\" \
             viewBox=\"0 0 {w} {h}\" font-family=\"Helvetica, Arial, sans-serif\">\n\
             <rect width=\"{w}\" height=\"{h}\" fill=\"white\"/>\n{}</svg>\n",
            self.body,
            w = n(width),
            h = n(height),
        )
    }

    /// `attributes` is written verbatim, such as `fill="none" stroke="black"`.
    pub fn polyline(&mut self, points: impl IntoIterator<Item = (f64, f64)>, attributes: &str) {
        let points: Vec<String> = points
            .into_iter()
            .map(|(x, y)| format!("{},{}", n(x), n(y)))
            .collect();
        let _ = writeln!(
            self.body,
            "<polyline points=\"{}\" {attributes}/>",
            points.join(" ")
        );
    }

    pub fn polygon(&mut self, points: impl IntoIterator<Item = (f64, f64)>, attributes: &str) {
        let points: Vec<String> = points
            .into_iter()
            .map(|(x, y)| format!("{},{}", n(x), n(y)))
            .collect();
        let _ = writeln!(
            self.body,
            "<polygon points=\"{}\" {attributes}/>",
            points.join(" ")
        );
    }

    pub fn line(&mut self, from: (f64, f64), to: (f64, f64), attributes: &str) {
        let _ = writeln!(
            self.body,
            "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" {attributes}/>",
            n(from.0),
            n(from.1),
            n(to.0),
            n(to.1)
        );
    }

    pub fn circle(&mut self, centre: (f64, f64), radius: f64, attributes: &str) {
        let _ = writeln!(
            self.body,
            "<circle cx=\"{}\" cy=\"{}\" r=\"{}\" {attributes}/>",
            n(centre.0),
            n(centre.1),
            n(radius)
        );
    }

    /// A path from raw path data, such as `M0,0 H10`.
    pub fn path(&mut self, data: &str, attributes: &str) {
        let _ = writeln!(self.body, "<path d=\"{data}\" {attributes}/>");
    }

    pub fn text(&mut self, at: (f64, f64), size: f64, text: &str, attributes: &str) {
        let _ = writeln!(
            self.body,
            "<text x=\"{}\" y=\"{}\" font-size=\"{}\" {attributes}>{}</text>",
            n(at.0),
            n(at.1),
            n(size),
            escape(text)
        );
    }

    /// Opens a group translated by `offset`; `end` closes it.
    pub fn group(&mut self, offset: (f64, f64), attributes: &str) {
        let _ = writeln!(
            self.body,
            "<g transform=\"translate({},{})\" {attributes}>",
            n(offset.0),
            n(offset.1)
        );
    }

    pub fn end(&mut self) {
        self.body.push_str("</g>\n");
    }

    /// A clip path named `id`: a circle about the origin.
    pub fn circle_clip(&mut self, id: &str, radius: f64) {
        let _ = writeln!(
            self.body,
            "<clipPath id=\"{id}\"><circle r=\"{}\"/></clipPath>",
            n(radius)
        );
    }
}

/// A coordinate to two decimals, without trailing zeros.
pub fn n(value: f64) -> String {
    let text = format!("{value:.2}");
    let text = text.trim_end_matches('0').trim_end_matches('.');
    if text == "-0" {
        "0".into()
    } else {
        text.into()
    }
}

fn escape(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}
