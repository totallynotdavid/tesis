//! The Fortran chains, fed with files the Rust writers produced, must give the
//! golden outputs byte for byte. Needs gfortran; `reference/regenerate-golden`
//! builds the programs.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use kk::formats::{structure, waveforms};

fn reference() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("reference")
}

#[test]
fn fortran_reads_rust_written_inputs_exactly() {
    let examples = reference().join("examples/input");
    let inputs = tempfile::tempdir().unwrap();
    for name in ["fort.1", "fort.1.Peru"] {
        let text = fs::read_to_string(examples.join(name)).unwrap();
        let parsed = waveforms::parse(&text).unwrap();
        fs::write(inputs.path().join(name), waveforms::write(&parsed)).unwrap();
    }
    for name in ["fort.2", "fort.2.Peru"] {
        let text = fs::read_to_string(examples.join(name)).unwrap();
        let parsed = structure::parse(&text).unwrap();
        fs::write(inputs.path().join(name), structure::write(&parsed)).unwrap();
    }

    let out = tempfile::tempdir().unwrap();
    let status = Command::new(reference().join("regenerate-golden"))
        .arg("--inputs")
        .arg(inputs.path())
        .arg("--out")
        .arg(out.path())
        .status()
        .unwrap();
    assert!(status.success());

    let golden = reference().join("golden");
    let mut compared = 0;
    for case in fs::read_dir(&golden).unwrap() {
        let case = case.unwrap().path();
        for file in fs::read_dir(&case).unwrap() {
            let file = file.unwrap().path();
            let ours = out
                .path()
                .join(case.file_name().unwrap())
                .join(file.file_name().unwrap());
            assert!(
                fs::read(&file).unwrap() == fs::read(&ours).unwrap(),
                "{} differs from the golden",
                ours.display()
            );
            compared += 1;
        }
    }
    assert_eq!(compared, 18);
}
