//! The figures `graphics` draws from `fort.4`, as SVG: the solution (source
//! time function, grid and mechanisms) and the waveform fits.

pub mod sphere;
mod svg;

use crate::formats::inversion::Shape;
use crate::formats::plot::{Plot, Subevent, Trace};
use crate::formats::solution::Mechanism;
use svg::Svg;

const WIDTH: f64 = 800.0;
const MARGIN: f64 = 30.0;
const BLACK: &str = "#111";
const SYNTHETIC: &str = "#d62728";

/// The moment rate, the source grid and one mechanism per subevent.
pub fn solution(p: &Plot) -> String {
    let mut svg = Svg::new();
    header(&mut svg, p);
    moment_rate(&mut svg, p, (MARGIN + 30.0, 100.0), (520.0, 110.0));
    if p.points.len() > 1 {
        grid(&mut svg, p, (680.0, 160.0), 60.0);
    }
    let height = mechanisms(&mut svg, p, 290.0);
    svg.finish(WIDTH, height)
}

/// Observed and synthetic records, eight stations to a column, each phase
/// type on a common amplitude scale as `graphics` draws with `norm` 2.
pub fn waveforms(p: &Plot) -> String {
    const ROWS: usize = 8;
    const LABEL: f64 = 110.0;
    const TRACE: f64 = 280.0;
    const ROW: f64 = 68.0;
    const TOP: f64 = 140.0;
    let mut svg = Svg::new();
    header(&mut svg, p);
    legend(&mut svg, (MARGIN, 92.0));
    let scales = phase_scales(&p.stations);
    let duration = p.source.len() as f64 * p.dt as f64;
    let columns = p.stations.len().div_ceil(ROWS).max(1);
    for (c, chunk) in p.stations.chunks(ROWS).enumerate() {
        let left = MARGIN + c as f64 * (LABEL + TRACE + 20.0) + LABEL;
        time_axis(&mut svg, (left, TOP - 14.0), TRACE, duration, true);
        for (r, st) in chunk.iter().enumerate() {
            let middle = TOP + 20.0 + r as f64 * ROW;
            trace(
                &mut svg,
                st,
                (left, middle),
                TRACE,
                scales[st.phase as usize],
            );
        }
    }
    let width = MARGIN * 2.0 + columns as f64 * (LABEL + TRACE + 20.0) - 20.0;
    let rows = p.stations.len().min(ROWS) as f64;
    svg.finish(width.max(WIDTH), TOP + rows * ROW + 10.0)
}

fn header(svg: &mut Svg, p: &Plot) {
    svg.text((MARGIN, 42.0), 22.0, &p.title, "font-weight=\"bold\"");
    let t = p.time_function;
    let line = format!(
        "{} time function, t1 {} s, t2 {} s · {} subevents · residual variance {:.4}",
        shape(t.shape),
        t.t1,
        t.t2,
        p.subevents.len(),
        p.variance
    );
    svg.text((MARGIN, 66.0), 13.0, &line, "fill=\"#444\"");
}

fn shape(s: Shape) -> &'static str {
    match s {
        Shape::Impulse => "impulse",
        Shape::Trapezoid => "trapezoid",
        Shape::Cosine => "cosine",
        Shape::Triangle => "triangle",
    }
}

/// Each subevent's time function, scaled by its moment and shifted to its
/// onset, shaded, and their sum outlined.
fn moment_rate(svg: &mut Svg, p: &Plot, origin: (f64, f64), size: (f64, f64)) {
    let (left, top) = origin;
    let (width, height) = size;
    let nt = p.source.len();
    let dt = p.dt as f64;
    let rates: Vec<Vec<f64>> = p
        .subevents
        .iter()
        .map(|e| {
            shifted(
                &p.source,
                (e.time as f64 / dt).round() as isize,
                e.mechanism,
            )
        })
        .collect();
    let total: Vec<f64> = (0..nt).map(|i| rates.iter().map(|r| r[i]).sum()).collect();
    let peak = total.iter().cloned().fold(0.0, f64::max);
    let peak = if peak > 0.0 { peak } else { 1.0 };
    let base = top + height;
    let x = |i: usize| left + width * i as f64 / nt as f64;
    let y = |v: f64| base - height * v / peak;
    svg.text(
        (left, top - 8.0),
        13.0,
        "Moment rate",
        "font-weight=\"bold\"",
    );
    for rate in &rates {
        let outline = (0..nt).map(|i| (x(i), y(rate[i])));
        let closed = outline.chain([(x(nt - 1), base), (x(0), base)]);
        svg.polygon(
            closed,
            "fill=\"#999\" fill-opacity=\"0.45\" stroke=\"none\"",
        );
    }
    svg.polyline(
        (0..nt).map(|i| (x(i), y(total[i]))),
        &format!("fill=\"none\" stroke=\"{BLACK}\" stroke-width=\"1.5\""),
    );
    time_axis(svg, (left, base), width, nt as f64 * dt, false);
    let moment = p.total.moment as f64;
    let mw = (moment.log10() + 18.0 - 9.1) / 1.5;
    let label = format!("Mo {} N·m   Mw {mw:.2}", scientific(moment * 1e18));
    svg.text(
        (left + width, top - 8.0),
        13.0,
        &label,
        "text-anchor=\"end\"",
    );
}

fn shifted(source: &[f32], shift: isize, m: Mechanism) -> Vec<f64> {
    (0..source.len() as isize)
        .map(|i| match usize::try_from(i - shift) {
            Ok(k) if k < source.len() => source[k] as f64 * m.moment as f64,
            _ => 0.0,
        })
        .collect()
}

/// A time axis of `duration` seconds, ticks below it, or above when `above`.
fn time_axis(svg: &mut Svg, origin: (f64, f64), width: f64, duration: f64, above: bool) {
    let (left, y) = origin;
    let stroke = format!("stroke=\"{BLACK}\" stroke-width=\"0.8\"");
    svg.line((left, y), (left + width, y), &stroke);
    let step = tick_step(duration);
    let (tick, text) = if above { (-4.0, -7.0) } else { (4.0, 16.0) };
    let mut t = 0.0;
    while t <= duration + 1e-6 {
        let x = left + width * t / duration;
        svg.line((x, y), (x, y + tick), &stroke);
        svg.text(
            (x, y + text),
            10.0,
            &format!("{t}"),
            "text-anchor=\"middle\"",
        );
        t += step;
    }
    svg.text((left + width + 6.0, y + text), 10.0, "s", "");
}

/// The smallest round step that gives at most six intervals.
fn tick_step(duration: f64) -> f64 {
    [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0]
        .into_iter()
        .find(|&s| duration / s <= 6.0)
        .unwrap_or(1000.0)
}

/// The searched points and the subevents in map view, north up, about the
/// hypocentre.
fn grid(svg: &mut Svg, p: &Plot, centre: (f64, f64), radius: f64) {
    let (cx, cy) = centre;
    let far = p
        .points
        .iter()
        .map(|&(r, _)| r)
        .chain(p.subevents.iter().map(|e| e.distance))
        .fold(0.0, f32::max) as f64;
    let scale = if far > 0.0 { radius / far } else { 1.0 };
    let at = |r: f32, azimuth: f32| {
        let (s, c) = (azimuth as f64).to_radians().sin_cos();
        (cx + r as f64 * scale * s, cy - r as f64 * scale * c)
    };
    svg.text(
        (cx - radius, cy - radius - 18.0),
        13.0,
        "Grid",
        "font-weight=\"bold\"",
    );
    svg.polyline(
        p.points.iter().map(|&(r, a)| at(r, a)),
        "fill=\"none\" stroke=\"#999\" stroke-width=\"0.8\"",
    );
    for &(r, a) in &p.points {
        svg.circle(at(r, a), 3.0, "fill=\"white\" stroke=\"#666\"");
    }
    let cross = format!("stroke=\"{BLACK}\" stroke-width=\"1\"");
    svg.line((cx - 6.0, cy), (cx + 6.0, cy), &cross);
    svg.line((cx, cy - 6.0), (cx, cy + 6.0), &cross);
    for (n, e) in p.subevents.iter().enumerate() {
        let (x, y) = at(e.distance, e.azimuth);
        svg.circle((x, y), 3.5, &format!("fill=\"{BLACK}\""));
        svg.text((x + 5.0, y - 5.0), 10.0, &(n + 1).to_string(), "");
    }
    let (nx, ny) = (cx + radius + 4.0, cy - radius + 8.0);
    svg.polygon(
        [(nx, ny - 12.0), (nx - 4.0, ny), (nx + 4.0, ny)],
        &format!("fill=\"{BLACK}\""),
    );
    svg.text((nx, ny + 12.0), 11.0, "N", "text-anchor=\"middle\"");
    let km = format!("{} km", (far.round()) as i64);
    svg.text(
        (cx, cy + radius + 18.0),
        10.0,
        &km,
        "text-anchor=\"middle\" fill=\"#444\"",
    );
    svg.line(
        (cx, cy + radius + 6.0),
        (cx + radius, cy + radius + 6.0),
        &cross,
    );
}

/// One focal sphere per subevent in onset order, sized by the square root of
/// its moment, then the total. Returns the height used.
fn mechanisms(svg: &mut Svg, p: &Plot, top: f64) -> f64 {
    const PER_ROW: usize = 5;
    const CELL: (f64, f64) = (148.0, 190.0);
    const LARGEST: f64 = 55.0;
    svg.text(
        (MARGIN, top - 10.0),
        13.0,
        "Mechanisms",
        "font-weight=\"bold\"",
    );
    let largest = p
        .subevents
        .iter()
        .map(|e| e.mechanism.moment.abs())
        .fold(0.0, f32::max);
    let mut order: Vec<usize> = (0..p.subevents.len()).collect();
    order.sort_by(|&a, &b| p.subevents[a].time.total_cmp(&p.subevents[b].time));
    let mut balls: Vec<Ball> = order
        .iter()
        .map(|&n| {
            let e = &p.subevents[n];
            let share = if largest > 0.0 {
                (e.mechanism.moment.abs() / largest).sqrt() as f64
            } else {
                1.0
            };
            Ball::subevent(n, e, (share * LARGEST).max(6.0))
        })
        .collect();
    if p.subevents.len() > 1 {
        balls.push(Ball {
            label: "Total".into(),
            note: None,
            mechanism: p.total,
            tensor: p.tensor,
            radius: LARGEST,
        });
    }
    for (k, ball) in balls.iter().enumerate() {
        let (row, column) = (k / PER_ROW, k % PER_ROW);
        let centre = (
            MARGIN + CELL.0 * (column as f64 + 0.5),
            top + 25.0 + CELL.1 * row as f64 + LARGEST,
        );
        ball.draw(svg, centre, &p.stations, &format!("ball{k}"));
    }
    top + CELL.1 * balls.len().div_ceil(PER_ROW) as f64 + 10.0
}

struct Ball {
    label: String,
    note: Option<String>,
    mechanism: Mechanism,
    tensor: [f32; 9],
    radius: f64,
}

impl Ball {
    fn subevent(n: usize, e: &Subevent, radius: f64) -> Self {
        Self {
            label: format!("{}  at {:.1} s", n + 1, e.time),
            note: Some(format!("{:.0} km, {:.0}°", e.distance, e.azimuth)),
            mechanism: e.mechanism,
            tensor: e.tensor,
            radius,
        }
    }

    fn draw(&self, svg: &mut Svg, centre: (f64, f64), stations: &[Trace], id: &str) {
        const PLANE_POINTS: usize = 90;
        let r = self.radius;
        let m = self.mechanism;
        svg.group(centre, "");
        svg.circle_clip(id, r);
        // Half-pixel rows keep the stairs of a shallow boundary under the
        // antialiasing.
        let rows = (4.0 * r).ceil() as usize;
        let mut fill = String::new();
        for (up, left, right) in sphere::compressions(&sphere::tensor(&self.tensor), rows) {
            let y = -(up + 1.0 / rows as f64) * r;
            fill.push_str(&format!(
                "M{},{}H{}v{}H{}Z",
                svg::n(left * r),
                svg::n(y),
                svg::n(right * r),
                svg::n(2.0 * r / rows as f64),
                svg::n(left * r)
            ));
        }
        if !fill.is_empty() {
            svg.path(&fill, &format!("fill=\"#555\" clip-path=\"url(#{id})\""));
        }
        let (normal, slip) = sphere::normal_and_slip(m);
        for pole in [normal, slip] {
            svg.polyline(
                sphere::great_circle(pole, PLANE_POINTS)
                    .into_iter()
                    .map(|(x, y)| (x * r, -y * r)),
                &format!("fill=\"none\" stroke=\"{BLACK}\" stroke-width=\"0.9\""),
            );
        }
        svg.circle(
            (0.0, 0.0),
            r,
            &format!("fill=\"none\" stroke=\"{BLACK}\" stroke-width=\"1.2\""),
        );
        for st in stations {
            let ray = sphere::ray(
                (st.azimuth as f64).to_radians(),
                (st.ray as f64).clamp(-1.0, 1.0).asin(),
            );
            let (x, y) = sphere::project(ray);
            svg.circle(
                (x * r, -y * r),
                2.2,
                &format!("fill=\"white\" stroke=\"{BLACK}\" stroke-width=\"0.7\""),
            );
        }
        let top = -r - 8.0;
        svg.text(
            (0.0, top),
            12.0,
            &self.label,
            "text-anchor=\"middle\" font-weight=\"bold\"",
        );
        let below = format!("{:.0}/{:.0}/{:.0}", m.strike, m.dip, m.slip);
        svg.text((0.0, r + 16.0), 11.0, &below, "text-anchor=\"middle\"");
        let moment = format!("Mo {} N·m", scientific(m.moment as f64 * 1e18));
        svg.text(
            (0.0, r + 30.0),
            10.0,
            &moment,
            "text-anchor=\"middle\" fill=\"#444\"",
        );
        if let Some(note) = &self.note {
            svg.text(
                (0.0, r + 43.0),
                10.0,
                note,
                "text-anchor=\"middle\" fill=\"#444\"",
            );
        }
        svg.end();
    }
}

/// The largest peak-to-peak record of each phase type, by `ib`.
fn phase_scales(stations: &[Trace]) -> [f64; 5] {
    let mut scales = [0.0; 5];
    for st in stations {
        let k = (st.phase as usize).min(4);
        scales[k] = f64::max(scales[k], peak_to_peak(&st.observed));
    }
    scales.map(|s| if s > 0.0 { s } else { 1.0 })
}

fn peak_to_peak(x: &[f32]) -> f64 {
    let max = x.iter().cloned().fold(f32::MIN, f32::max);
    let min = x.iter().cloned().fold(f32::MAX, f32::min);
    if x.is_empty() {
        0.0
    } else {
        (max - min) as f64
    }
}

fn trace(svg: &mut Svg, st: &Trace, origin: (f64, f64), width: f64, scale: f64) {
    const HEIGHT: f64 = 50.0;
    let (left, middle) = origin;
    let nt = st.observed.len().max(1);
    let x = |i: usize| left + width * i as f64 / nt as f64;
    let y = |v: f32| middle - HEIGHT * v as f64 / scale;
    let phase = ["?", "P", "SV", "SH", "PP"][(st.phase as usize).min(4)];
    let azimuth = if st.azimuth < 0.0 {
        st.azimuth + 360.0
    } else {
        st.azimuth
    };
    let label = left - 104.0;
    svg.text(
        (label, middle - 4.0),
        13.0,
        st.name.trim(),
        "font-weight=\"bold\"",
    );
    let line = format!("{phase}  {azimuth:.1}°");
    svg.text((label, middle + 10.0), 10.0, &line, "fill=\"#444\"");
    let amplitude = scientific(peak_to_peak(&st.observed));
    svg.text((label, middle + 22.0), 10.0, &amplitude, "fill=\"#444\"");
    svg.polyline(
        st.observed.iter().enumerate().map(|(i, &v)| (x(i), y(v))),
        &format!("fill=\"none\" stroke=\"{BLACK}\" stroke-width=\"1.3\""),
    );
    svg.polyline(
        st.observed
            .iter()
            .zip(&st.residual)
            .enumerate()
            .map(|(i, (&o, &r))| (x(i), y(o - r))),
        &format!("fill=\"none\" stroke=\"{SYNTHETIC}\" stroke-width=\"1\""),
    );
}

fn legend(svg: &mut Svg, origin: (f64, f64)) {
    let (x, y) = origin;
    let observed = format!("stroke=\"{BLACK}\" stroke-width=\"1.3\"");
    let synthetic = format!("stroke=\"{SYNTHETIC}\" stroke-width=\"1\"");
    svg.line((x, y - 4.0), (x + 24.0, y - 4.0), &observed);
    svg.text((x + 30.0, y), 11.0, "observed", "");
    svg.line((x + 100.0, y - 4.0), (x + 124.0, y - 4.0), &synthetic);
    svg.text((x + 130.0, y), 11.0, "synthetic", "");
    let note = "peak-to-peak amplitude under each name; one scale per phase type";
    svg.text((x + 210.0, y), 11.0, note, "fill=\"#444\"");
}

/// Three significant figures, as `1.23e20`.
fn scientific(value: f64) -> String {
    format!("{value:.2e}")
}
