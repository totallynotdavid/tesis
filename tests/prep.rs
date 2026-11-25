//! `prep` against `conv.sac.farm`. The test builds a copy of the reference
//! program whose `fort.22` formats print every computed real to nine
//! significant digits, which identifies an `f32`. Both sides convert the
//! same synthetic SAC records, and every value must match exactly. Needs
//! gfortran.

use std::collections::VecDeque;
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::str::FromStr;

use kk::fft::Complex;
use kk::formats::waveforms::{Instrument, Station};
use kk::prep::earth::{Earth, Phase, Tables};
use kk::prep::{self, Event, Quantity, Record, Response, Skip, Time, Window};

/// The 2013-09-25 Yauca origin, 16:42:43.5 UTC on day 268.
const EVENT: (f32, f32) = (-15.84, -74.51);
const ORIGIN_MINUTE: i64 = 16 * 60 + 42;
const ORIGIN_SECOND: f32 = 43.5;
const DT0: f32 = 0.05;
const SAMPLES: usize = 4000;
/// Records start this many seconds before the phase onset.
const PRE_ONSET: f64 = 60.0;

struct Site {
    code: &'static str,
    channel: &'static str,
    latitude: f32,
    longitude: f32,
    phase: Phase,
    low: f32,
    high: f32,
    duration: f32,
    units: Quantity,
}

/// IU and II stations. TRQA is closer than 30 degrees and is skipped; PAB's
/// window runs past the end of its record.
const SITES: [Site; 7] = [
    site(
        "ANMO",
        "BHZ",
        34.95,
        -106.46,
        Phase::P,
        0.01,
        1.0,
        60.0,
        Quantity::Velocity,
    ),
    site(
        "HRV",
        "BHN",
        42.51,
        -71.56,
        Phase::P,
        0.02,
        0.5,
        60.0,
        Quantity::Displacement,
    ),
    site(
        "PAB",
        "BHE",
        39.55,
        -4.35,
        Phase::P,
        0.01,
        1.0,
        200.0,
        Quantity::Acceleration,
    ),
    site(
        "SJG",
        "BH1",
        18.11,
        -66.15,
        Phase::S,
        0.01,
        0.2,
        80.0,
        Quantity::Velocity,
    ),
    site(
        "RCBR",
        "BH2",
        -5.83,
        -35.90,
        Phase::P,
        0.0,
        1.0,
        60.0,
        Quantity::Displacement,
    ),
    site(
        "PFO",
        "BHZ",
        33.61,
        -116.46,
        Phase::P,
        0.01,
        0.0,
        60.0,
        Quantity::Acceleration,
    ),
    site(
        "TRQA",
        "BHZ",
        -38.06,
        -61.98,
        Phase::P,
        0.01,
        1.0,
        60.0,
        Quantity::Velocity,
    ),
];

#[allow(clippy::too_many_arguments)]
const fn site(
    code: &'static str,
    channel: &'static str,
    latitude: f32,
    longitude: f32,
    phase: Phase,
    low: f32,
    high: f32,
    duration: f32,
    units: Quantity,
) -> Site {
    Site {
        code,
        channel,
        latitude,
        longitude,
        phase,
        low,
        high,
        duration,
        units,
    }
}

impl Site {
    fn file(&self) -> String {
        format!("{}.{}.sac", self.code, self.channel)
    }

    fn window(&self) -> Window {
        let (lead, dt) = match self.phase {
            Phase::P => (10.0, 0.1),
            Phase::S => (20.0, 0.2),
        };
        Window {
            phase: self.phase,
            lead,
            duration: self.duration,
            dt,
            low: self.low,
            high: self.high,
        }
    }
}

fn reference() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("reference")
}

#[test]
fn prep_matches_conv_sac_farm() {
    let dir = tempfile::tempdir().unwrap();
    let conv = build_oracle(dir.path());
    let input = reference().join("examples/input");
    let table = |name: &str| fs::read_to_string(input.join(name)).unwrap();
    let tables = Tables::read(&table("jb.ptime"), &table("jb.stime"), &table("jb.table")).unwrap();
    for name in ["jb.ptime", "jb.stime", "jb.table"] {
        fs::copy(input.join(name), dir.path().join(name)).unwrap();
    }

    let responses: String = SITES.iter().map(response_entry).collect();
    fs::write(dir.path().join("response"), responses).unwrap();

    let events = [
        (Quantity::Displacement, 40.0),
        (Quantity::Velocity, 150.0),
        (Quantity::Acceleration, 600.0),
    ];
    for (output, depth) in events {
        let earth = tables.at(depth).unwrap();
        let event = Event {
            title: "2013 Yauca test".to_owned(),
            latitude: EVENT.0,
            longitude: EVENT.1,
            origin: Time {
                minute: ORIGIN_MINUTE,
                second: ORIGIN_SECOND,
            },
            output,
            distances: 30.0..=90.0,
        };
        let mut stdin = format!(
            "{}\n{} {} {depth} 16 42 {ORIGIN_SECOND} {} 30 90\n",
            event.title, EVENT.0, EVENT.1, output as i32
        );
        let mut ours = Vec::new();
        for (n, site) in SITES.iter().enumerate() {
            let record = record(site, &event, &earth, n as u64);
            fs::write(dir.path().join(site.file()), sac(site, &record)).unwrap();
            let window = site.window();
            let _ = writeln!(
                stdin,
                "{}\n{}\n{} {} {} {} {} {} 1 {}",
                site.file(),
                site.code,
                window.lead,
                window.duration,
                window.dt,
                window.low,
                window.high,
                window.phase as i32 + 1,
                site.units as i32
            );
            match prep::convert(&event, &earth, &record, &window) {
                Ok(station) => ours.push(station),
                Err(skip) => {
                    assert_eq!(site.code, "TRQA", "{}: {skip}", site.code);
                    assert!(matches!(skip, Skip::Distance(_)));
                }
            }
        }
        stdin.push_str("dummy\ndummy\n");
        assert_eq!(ours.len(), SITES.len() - 1);

        let fort22 = run(&conv, dir.path(), &stdin);
        let mut oracle = Oracle::new(&fort22);
        assert_eq!(oracle.line(), prep::title(&event, &earth));
        for station in &ours {
            let theirs = oracle.station(output);
            compare(station, &theirs, depth);
        }
        assert!(oracle.lines.next().is_none(), "fort.22 has more stations");
    }
}

/// Compiles `conv.sac.farm` with full-precision output formats.
fn build_oracle(dir: &Path) -> PathBuf {
    let source = fs::read_to_string(reference().join("main/conv.sac.farm.f")).unwrap();
    let replacements = [
        ("FORMAT(3F8.2,F8.3,F8.2,1x,i4)", "FORMAT(5es16.8,1x,i4)", 1),
        (
            "format(a12,i4,f7.2,f8.2,3f7.2)",
            "format(a12,i4,5es16.8)",
            1,
        ),
        ("'(e12.3,i7,f7.3)'", "'(es16.8,i7,es16.8)'", 1),
        ("'(10f7.3)'", "'(5es16.8)'", 2),
    ];
    let mut variant = source;
    for (from, to, count) in replacements {
        assert_eq!(variant.matches(from).count(), count, "{from}");
        variant = variant.replace(from, to);
    }
    assert!(variant.lines().all(|l| l.len() <= 72));
    let path = dir.join("conv.f");
    fs::write(&path, variant).unwrap();

    let makefile = fs::read_to_string(reference().join("Makefile")).unwrap();
    let flags = makefile
        .lines()
        .find_map(|l| l.strip_prefix("FFLAGS"))
        .and_then(|l| l.split_once("?="))
        .unwrap()
        .1;
    let conv = dir.join("conv");
    let status = Command::new("gfortran")
        .args(flags.split_whitespace())
        .arg(&path)
        .arg(reference().join("subroutine/sub.cfft.f"))
        .arg("-o")
        .arg(&conv)
        .status()
        .unwrap();
    assert!(status.success());
    conv
}

fn run(conv: &Path, dir: &Path, stdin: &str) -> String {
    let _ = fs::remove_file(dir.join("fort.22"));
    let mut child = Command::new(conv)
        .current_dir(dir)
        .stdin(Stdio::piped())
        .stdout(Stdio::null())
        .spawn()
        .unwrap();
    std::io::Write::write_all(&mut child.stdin.take().unwrap(), stdin.as_bytes()).unwrap();
    assert!(child.wait().unwrap().success());
    fs::read_to_string(dir.join("fort.22")).unwrap()
}

/// A record that starts [`PRE_ONSET`] seconds, give or take, before the
/// phase onset: seeded noise, then a decaying pulse at the onset.
fn record(site: &Site, event: &Event, earth: &Earth, seed: u64) -> Record {
    let geometry = earth.geometry(
        (event.latitude, event.longitude),
        (site.latitude, site.longitude),
    );
    let onset = f64::from(ORIGIN_SECOND) + f64::from(earth.travel(site.phase, geometry.distance));
    let start = (onset - PRE_ONSET).floor() as i64;
    let mut state = seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) | 1;
    let samples = (0..SAMPLES)
        .map(|i| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let noise = ((state >> 40) as f32 / 16_777_216.0 - 0.5) * 200.0;
            let t = i as f64 * f64::from(DT0) - (PRE_ONSET - 0.25);
            let pulse = if t > 0.0 {
                5000.0 * (t * std::f64::consts::TAU / 8.0).sin() * (-t / 10.0).exp()
            } else {
                0.0
            };
            noise + pulse as f32
        })
        .collect();
    Record {
        code: site.code.to_owned(),
        channel: site.channel.to_owned(),
        latitude: site.latitude,
        longitude: site.longitude,
        start: Time {
            minute: ORIGIN_MINUTE + start.div_euclid(60),
            second: start.rem_euclid(60) as f32 + 250.0 / 1000.0,
        },
        dt: DT0,
        samples,
        response: response(site.units),
    }
}

/// A broadband response; one zero at the origin per integration from the
/// sensor's units to displacement.
fn response(units: Quantity) -> Response {
    let mut zeros = vec![Complex::new(-31.63, 0.0)];
    zeros.resize(4 - units as usize, Complex::default());
    Response {
        units,
        normalisation: 6.0077e7,
        sensitivity: 2.5e9,
        poles: [
            (-0.03701, 0.03701),
            (-0.03701, -0.03701),
            (-251.3, 0.0),
            (-131.0, 467.3),
            (-131.0, -467.3),
        ]
        .map(|(re, im)| Complex::new(re, im))
        .to_vec(),
        zeros,
    }
}

/// One entry of the `response` file `conv.sac.farm` reads.
fn response_entry(site: &Site) -> String {
    let r = response(site.units);
    let mut out = format!(
        "{}\n{:<28}{}\n{:<28}{}\n{:<28}{}\n",
        site.file(),
        " A0 normalization factor:",
        r.normalisation,
        " Number of zeroes:",
        r.zeros.len(),
        " Number of poles:",
        r.poles.len()
    );
    for (i, z) in r.zeros.iter().chain(&r.poles).enumerate() {
        let _ = writeln!(out, "{i} {} {} 0 0", z.re, z.im);
    }
    let _ = writeln!(out, "{:<28}{}", " Sensitivity:", r.sensitivity);
    out
}

/// A SAC alphanumeric file holding the header fields `conv.sac.farm` reads.
fn sac(site: &Site, record: &Record) -> String {
    let mut floats = [-12345.0f32; 70];
    floats[0] = record.dt;
    floats[31] = record.latitude;
    floats[32] = record.longitude;
    floats[33] = 0.0;
    floats[34] = 0.0;
    let mut out = String::new();
    for row in floats.chunks(5) {
        for v in row {
            // A G15.7 field without a decimal point would be scaled by 1e-7.
            let s = v.to_string();
            let s = if s.contains('.') { s } else { s + "." };
            let _ = write!(out, "{s:>15}");
        }
        out.push('\n');
    }
    let second = record.start.second.trunc() as i64;
    let msec = ((record.start.second - second as f32) * 1000.0).round() as i64;
    let mut ints = [-12345i64; 40];
    ints[..6].copy_from_slice(&[
        2013,
        268,
        record.start.minute / 60,
        record.start.minute % 60,
        second,
        msec,
    ]);
    ints[6] = 6;
    ints[9] = record.samples.len() as i64;
    for row in ints.chunks(5) {
        for v in row {
            let _ = write!(out, "{v:>10}");
        }
        out.push('\n');
    }
    let _ = writeln!(out, "{:<8}{:<16}", site.code, "-12345");
    for _ in 0..5 {
        let _ = writeln!(out, "{:<8}{:<8}{:<8}", "-12345", "-12345", "-12345");
    }
    let _ = writeln!(out, "{:<8}{:<8}{:<8}", "-12345", "-12345", site.channel);
    let _ = writeln!(out, "{:<8}{:<8}{:<8}", "-12345", "-12345", "-12345");
    for row in record.samples.chunks(5) {
        let row: Vec<String> = row.iter().map(f32::to_string).collect();
        let _ = writeln!(out, "{}", row.join(" "));
    }
    out
}

/// `fort.22` from the variant, read back into stations.
struct Oracle<'a> {
    lines: std::str::Lines<'a>,
    tokens: VecDeque<&'a str>,
}

impl<'a> Oracle<'a> {
    fn new(text: &'a str) -> Self {
        Self {
            lines: text.lines(),
            tokens: VecDeque::new(),
        }
    }

    fn line(&mut self) -> &'a str {
        assert!(self.tokens.is_empty());
        self.lines.next().unwrap()
    }

    fn next<T: FromStr>(&mut self) -> T
    where
        T::Err: std::fmt::Debug,
    {
        while self.tokens.is_empty() {
            self.tokens
                .extend(self.lines.next().unwrap().split_whitespace());
        }
        self.tokens.pop_front().unwrap().parse().unwrap()
    }

    fn station(&mut self, output: Quantity) -> Station {
        let line = self.line();
        let (code, rest) = line.split_at(12);
        self.tokens.extend(rest.split_whitespace());
        let minute: i64 = self.next();
        let [second, lead, low, high] = [(); 4].map(|_| self.next::<f32>());
        let name = format!("{code:<12.12}{minute:>4}{second:>7.2}{lead:>8.2}{low:>7.2}{high:>7.2}");
        let [azimuth, back_azimuth, distance, slowness, spreading] =
            [(); 5].map(|_| self.next::<f32>());
        let onset = self.next();
        assert_eq!(self.next::<i32>(), 1);
        let (phase, component) = (self.next(), self.next());
        assert_eq!([self.next::<i32>(), self.next()], [0, 0]);
        // The reference writes power 1 for velocity and 0 otherwise; green
        // reads it as the number of differentiations, so prep writes 2 for
        // acceleration.
        assert_eq!(self.next::<i32>(), i32::from(output == Quantity::Velocity));
        assert_eq!([self.next::<f32>(), self.next()], [1.0, 1.0]);
        let (peak, count, dt): (f32, usize, f32) = (self.next(), self.next(), self.next());
        let samples = (0..count).map(|_| self.next()).collect();
        Station {
            name,
            azimuth,
            back_azimuth,
            distance,
            slowness,
            spreading,
            onset,
            phase,
            component,
            instrument: Instrument::Gdsn {
                power: output as i32,
                magnification: 1.0,
                normalisation: 1.0,
                poles: Vec::new(),
                zeros: Vec::new(),
            },
            peak,
            dt,
            samples,
        }
    }
}

fn compare(ours: &Station, theirs: &Station, depth: f32) {
    let at = format!("{} at {depth} km", ours.name);
    assert_eq!(ours.samples.len(), theirs.samples.len(), "{at}");
    let first = ours
        .samples
        .iter()
        .zip(&theirs.samples)
        .position(|(a, b)| a != b);
    if let Some(i) = first {
        panic!(
            "{at}: sample {i} is {} against {}",
            ours.samples[i], theirs.samples[i]
        );
    }
    let header = |s: &Station| Station {
        samples: Vec::new(),
        ..s.clone()
    };
    assert_eq!(header(ours), header(theirs), "{at}");
}
