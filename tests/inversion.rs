//! `inversion` on the Denali example must write the golden `fort.16`,
//! `fort.4` and standard output byte for byte, and must match the reference
//! program where the golden does not reach.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use kk::formats::inversion::{Basis, Grid, Parameters, Shape, TimeFunction};
use kk::formats::{greens, inversion, solution, structure, waveforms};
use kk::green::{self, Line, Sampling};

fn reference() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("reference")
}

fn read(path: &str) -> Vec<u8> {
    fs::read(reference().join(path)).unwrap()
}

fn text(path: &str) -> String {
    String::from_utf8(read(path)).unwrap()
}

#[test]
fn inversion_matches_the_denali_golden_files() {
    let waveforms = waveforms::parse(&text("examples/input/fort.1")).unwrap();
    let greens = greens::read(&read("golden/denali/fort.3")).unwrap();
    let parameters = inversion::parse(&text("examples/input/i_inversion")).unwrap();
    let ours = kk::inversion::invert(&parameters, &waveforms, &greens).unwrap();
    // The golden stdout holds every program's output, each after a header.
    let stdout = text("golden/denali/stdout");
    let section = stdout
        .split("== ")
        .find_map(|s| s.strip_prefix("inversion\n"))
        .unwrap();
    assert_eq!(solution::summary(&ours), section);
    assert_eq!(solution::report(&ours), text("golden/denali/fort.16"));
    assert!(
        solution::plot(&ours) == read("golden/denali/fort.4"),
        "same report, different fort.4"
    );
}

/// The golden runs one basis, one time function and one depth. This runs the
/// reference `inversion` on Denali under every basis, time function and grid
/// kind, with one and three depths, one to thirteen subevents, and stations
/// left out by weight. Needs gfortran.
#[test]
fn inversion_matches_the_reference_program_on_every_basis() {
    let dir = tempfile::tempdir().unwrap();
    let program = build_inversion(dir.path());
    let fort1 = text("examples/input/fort.1");
    let waveforms = waveforms::parse(&fort1).unwrap();
    let one_depth = read("golden/denali/fort.3");
    let three_depths = greens::write(&three_depths(&waveforms));
    let denali = inversion::parse(&text("examples/input/i_inversion")).unwrap();
    let points = Grid::Points(vec![(0.0, 0.0), (30.0, 120.0), (60.0, 300.0)]);
    let line = Grid::Line {
        azimuth: 120.0,
        spacing: 20.0,
        count: 4,
        reference: 2,
    };
    let shape = |shape, t1, t2| TimeFunction { shape, t1, t2 };
    let cases = [
        (
            Basis::Deviatoric,
            shape(Shape::Impulse, 0.0, 0.0),
            1,
            &points,
            &one_depth,
        ),
        (
            Basis::DoubleCouple,
            shape(Shape::Trapezoid, 6.0, 4.0),
            2,
            &line,
            &three_depths,
        ),
        (
            Basis::StrikeSlip,
            shape(Shape::Cosine, 4.0, 8.0),
            3,
            &points,
            &three_depths,
        ),
        (
            Basis::MomentTensor,
            shape(Shape::Cosine, 6.0, 3.0),
            2,
            &line,
            &one_depth,
        ),
        (
            Basis::VerticalDoubleCouple,
            shape(Shape::Triangle, 3.0, 8.0),
            13,
            &points,
            &three_depths,
        ),
        (
            Basis::MomentTensor,
            shape(Shape::Trapezoid, 2.0, 5.0),
            1,
            &line,
            &three_depths,
        ),
    ];
    for (n, (basis, time_function, subevents, grid, fort3)) in cases.into_iter().enumerate() {
        let mut parameters = Parameters {
            length: 60.0,
            time_function,
            spans: (0..subevents).map(|e| 20.0 + 5.0 * e as f32).collect(),
            basis,
            grid: grid.clone(),
            ..denali.clone()
        };
        if n == 5 {
            // Ten weights for the kept stations, two of them zero.
            parameters.weights = vec![1.0, 0.0, 1.0, 0.5, 1.0, 1.0, 0.0, 1.0, 2.0, 1.0];
            parameters.shift = 0;
        }
        let case = format!("case {n}: {basis:?}, {time_function:?}, {subevents} subevents");
        fs::write(dir.path().join("fort.1"), &fort1).unwrap();
        fs::write(dir.path().join("fort.3"), fort3).unwrap();
        fs::write(
            dir.path().join("i_inversion"),
            inversion::write(&parameters),
        )
        .unwrap();
        let output = Command::new(&program)
            .current_dir(dir.path())
            .output()
            .unwrap();
        assert!(output.status.success(), "{case}");
        let theirs = |name: &str| fs::read(dir.path().join(name)).unwrap();
        let greens = greens::read(fort3).unwrap();
        let ours = kk::inversion::invert(&parameters, &waveforms, &greens).unwrap();
        assert_eq!(
            solution::summary(&ours),
            String::from_utf8(output.stdout).unwrap(),
            "{case}"
        );
        assert_eq!(
            solution::report(&ours),
            String::from_utf8(theirs("fort.16")).unwrap(),
            "{case}"
        );
        assert!(solution::plot(&ours) == theirs("fort.4"), "{case}: fort.4");
    }
}

/// Denali's base elements at three depths, 10 km apart about 15 km.
fn three_depths(waveforms: &waveforms::Waveforms) -> greens::Greens {
    let structure = structure::parse(&text("examples/input/fort.2")).unwrap();
    let sampling = Sampling { nt: 256, dt: 0.5 };
    let line = Line {
        h0: 15.0,
        dk: 10.0,
        nk: 3,
        k0: 2,
        dip: 90.0,
    };
    green::base(&structure, waveforms, sampling, &line).unwrap()
}

/// The reference `inversion`, built with the Makefile's flags from its
/// sources.
fn build_inversion(dir: &Path) -> PathBuf {
    let makefile = text("Makefile");
    let flags = makefile
        .lines()
        .find_map(|l| l.strip_prefix("FFLAGS"))
        .and_then(|l| l.split_once("?="))
        .unwrap()
        .1;
    let subroutines = [
        "conv", "corr", "source", "zmom", "det0", "det00", "det1", "mtrx", "cubicd", "quadrd",
        "dcp", "eig1",
    ];
    let recipes = ["jacobi", "eigsrt", "ludcmp", "lubksb"];
    let program = dir.join("inversion");
    let status = Command::new("gfortran")
        .args(flags.split_whitespace())
        .arg(reference().join("main/inversion.f"))
        .args(subroutines.map(|s| reference().join(format!("subroutine/sub.{s}.f"))))
        .args(recipes.map(|s| reference().join(format!("recipes/{s}.f"))))
        .arg("-o")
        .arg(&program)
        .status()
        .unwrap();
    assert!(status.success());
    program
}
