//! `green` on the Denali example and `green3` on the Peru example must write
//! the golden `fort.3` and `fort.33` byte for byte, and `green` must match
//! the reference program where the goldens do not reach.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use kk::formats::waveforms::Instrument;
use kk::formats::{fixed, greens, structure, waveforms};
use kk::green::{self, Line, Mechanism, Sampling};

fn reference() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("reference")
}

fn inputs(fort1: &str, fort2: &str) -> (structure::Structure, waveforms::Waveforms) {
    let input = reference().join("examples/input");
    let read = |name: &str| fs::read_to_string(input.join(name)).unwrap();
    (
        structure::parse(&read(fort2)).unwrap(),
        waveforms::parse(&read(fort1)).unwrap(),
    )
}

#[test]
fn green_matches_the_denali_golden_file() {
    let (structure, waveforms) = inputs("fort.1", "fort.2");
    // `i_green`: 256 .5 15 5 1 1 90
    let sampling = Sampling { nt: 256, dt: 0.5 };
    let line = Line {
        h0: 15.0,
        dk: 5.0,
        nk: 1,
        k0: 1,
        dip: 90.0,
    };
    let ours = green::base(&structure, &waveforms, sampling, &line).unwrap();
    let golden = fs::read(reference().join("golden/denali/fort.3")).unwrap();
    let theirs = greens::read(&golden).unwrap();
    let traces = |g: &greens::Greens| -> Vec<_> {
        g.stations
            .iter()
            .map(|s| (s.name.clone(), s.depths.iter().flatten().cloned().collect()))
            .collect()
    };
    assert_eq!(first_difference(&traces(&ours), &traces(&theirs)), None);
    assert!(
        greens::write(&ours) == golden,
        "same traces, different bytes"
    );
}

#[test]
fn green3_matches_the_peru_golden_file() {
    let (structure, waveforms) = inputs("fort.1.Peru", "fort.2.Peru");
    // `i_green3`: 256 0.5 309 21 59 30 30 6 2. The points run down the fault.
    let sampling = Sampling { nt: 256, dt: 0.5 };
    let mechanism = Mechanism::new(309.0, 21.0, 59.0);
    let line = Line {
        h0: 30.0,
        dk: 30.0,
        nk: 6,
        k0: 2,
        dip: mechanism.dip,
    };
    let ours = green::fixed(&structure, &waveforms, sampling, &line, mechanism).unwrap();
    let golden = fs::read(reference().join("golden/peru/fort.33")).unwrap();
    let theirs = fixed::read(&golden).unwrap();
    let traces = |g: &fixed::Greens| -> Vec<_> {
        g.stations
            .iter()
            .map(|s| (s.name.clone(), s.depths.clone()))
            .collect()
    };
    assert_eq!(first_difference(&traces(&ours), &traces(&theirs)), None);
    assert!(
        fixed::write(&ours) == golden,
        "same traces, different bytes"
    );
}

/// The goldens cover only P vertical and SH on GDSN records without poles.
/// This runs the reference `green` on Denali stations rewritten to every
/// other phase, component and instrument, over depths in each source layer,
/// on its boundary and in the half-space, under the Denali structure and the
/// Peru one with a water layer. Needs gfortran.
///
/// The reference `bodyw` computes `tp` for P and PP only, so its SV
/// synthetic uses the `tp` the previous call saved. `green` gives SV its
/// own, and the oracle's `bodyw` is patched to do the same.
#[test]
fn green_matches_the_reference_program_on_every_phase_and_instrument() {
    let dir = tempfile::tempdir().unwrap();
    let program = build_green(dir.path());
    let fort1 = waveforms::write(&variant_stations());
    let waveforms = waveforms::parse(&fort1).unwrap();
    let sampling = Sampling { nt: 256, dt: 0.5 };
    let line = Line {
        h0: 15.0,
        dk: 10.0,
        nk: 4,
        k0: 2,
        dip: 90.0,
    };
    let input = reference().join("examples/input");
    for fort2 in ["fort.2", "fort.2.Peru"] {
        let text = fs::read_to_string(input.join(fort2)).unwrap();
        let structure = structure::parse(&text).unwrap();
        fs::write(dir.path().join("fort.1"), &fort1).unwrap();
        fs::write(dir.path().join("fort.2"), &text).unwrap();
        fs::write(dir.path().join("i_green"), "256 .5 15 10 4 2 90\n").unwrap();
        let _ = fs::remove_file(dir.path().join("fort.3"));
        let status = Command::new(&program)
            .current_dir(dir.path())
            .stdout(Stdio::null())
            .status()
            .unwrap();
        assert!(status.success());
        let theirs = greens::read(&fs::read(dir.path().join("fort.3")).unwrap()).unwrap();
        let ours = green::base(&structure, &waveforms, sampling, &line).unwrap();
        let traces = |g: &greens::Greens| -> Vec<_> {
            g.stations
                .iter()
                .map(|s| (s.name.clone(), s.depths.iter().flatten().cloned().collect()))
                .collect()
        };
        assert_eq!(theirs.stations.len(), waveforms.stations.len());
        assert_eq!(
            first_difference(&traces(&ours), &traces(&theirs)),
            None,
            "{fort2}"
        );
    }
}

/// The first Denali stations, each given one (phase, component, instrument)
/// the goldens do not reach.
fn variant_stations() -> waveforms::Waveforms {
    let text = fs::read_to_string(reference().join("examples/input/fort.1")).unwrap();
    let mut waveforms = waveforms::parse(&text).unwrap();
    let wwssn = Instrument::Wwssn {
        magnification: 1500.0,
        seismometer_period: 15.0,
        galvanometer_period: 100.0,
        seismometer_damping: 0.93,
        galvanometer_damping: 1.0,
        coupling: 0.003,
    };
    let sts1 = Instrument::Gdsn {
        power: 0,
        magnification: 2.3e3,
        normalisation: 3948.58,
        poles: vec![
            (-0.01234, 0.01234),
            (-0.01234, -0.01234),
            (-39.18, 49.12),
            (-39.18, -49.12),
        ],
        zeros: vec![(0.0, 0.0), (0.0, 0.0)],
    };
    let velocity = Instrument::Gdsn {
        power: 1,
        magnification: 1.0,
        normalisation: 1.0,
        poles: Vec::new(),
        zeros: Vec::new(),
    };
    let cases = [
        (1, 2, &sts1),
        (1, 3, &wwssn),
        (2, 1, &sts1),
        (2, 2, &wwssn),
        (3, 1, &wwssn),
        (3, 1, &sts1),
        (4, 1, &velocity),
        (4, 2, &sts1),
        (4, 3, &wwssn),
        (1, 1, &wwssn),
        (2, 2, &velocity),
    ];
    waveforms.stations.truncate(cases.len());
    assert_eq!(waveforms.stations.len(), cases.len());
    for (station, (phase, component, instrument)) in waveforms.stations.iter_mut().zip(cases) {
        station.phase = phase;
        station.component = component;
        station.instrument = instrument.clone();
    }
    waveforms
}

/// The reference `green`, built with the Makefile's flags from its sources,
/// with SV given its own `tp` in `bodyw`.
fn build_green(dir: &Path) -> PathBuf {
    let makefile = fs::read_to_string(reference().join("Makefile")).unwrap();
    let flags = makefile
        .lines()
        .find_map(|l| l.strip_prefix("FFLAGS"))
        .and_then(|l| l.split_once("?="))
        .unwrap()
        .1;
    let source = reference().join("subroutine/sub.bodyw.f");
    let bodyw = fs::read_to_string(source).unwrap();
    let from = "          TS=YS*DH\n";
    assert_eq!(bodyw.matches(from).count(), 1);
    let own_tp = "          YP=SQRT(1/VP(LL)**2-P**2)\n          TP=YP*DH\n";
    let bodyw_path = dir.join("bodyw.f");
    fs::write(&bodyw_path, bodyw.replace(from, &format!("{from}{own_tp}"))).unwrap();
    let subroutines = [
        "refl", "cnvr", "reflsh", "cnvrsh", "radp", "qfm", "instg", "instw", "prod", "clear",
        "cfft",
    ];
    let program = dir.join("green");
    let status = Command::new("gfortran")
        .args(flags.split_whitespace())
        .arg(reference().join("main/green.f"))
        .arg(&bodyw_path)
        .args(subroutines.map(|s| reference().join(format!("subroutine/sub.{s}.f"))))
        .arg("-o")
        .arg(&program)
        .status()
        .unwrap();
    assert!(status.success());
    program
}

type Traces = Vec<(String, Vec<Vec<f32>>)>;

/// Where two sets of station traces first differ, bit for bit.
fn first_difference(a: &Traces, b: &Traces) -> Option<String> {
    if a.len() != b.len() {
        return Some(format!("{} stations != {}", a.len(), b.len()));
    }
    for ((na, ta), (nb, tb)) in a.iter().zip(b) {
        if na.trim_end() != nb.trim_end() {
            return Some(format!("station {na:?} != {nb:?}"));
        }
        for (t, (xa, xb)) in ta.iter().zip(tb).enumerate() {
            for (i, (x, y)) in xa.iter().zip(xb).enumerate() {
                if x.to_bits() != y.to_bits() {
                    return Some(format!("{na} trace {t} sample {i}: {x:e} != {y:e}"));
                }
            }
        }
    }
    None
}
