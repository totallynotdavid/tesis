//! The Kikuchi-Kanamori text files. Readers follow Fortran list-directed
//! rules; writers emit values that the Fortran programs read back exactly.

pub mod records;
pub mod structure;
pub mod waveforms;

use std::fmt::Write as _;

/// Reals right-aligned in 8 columns. `f32` Display prints the shortest
/// decimal that parses back to the same bits, and gfortran's list-directed
/// read rounds correctly, so the Fortran sees the exact value.
fn reals(values: &[f32]) -> String {
    let mut out = String::new();
    for value in values {
        let _ = write!(out, " {value:>7}");
    }
    out
}

fn write_rows(out: &mut String, values: &[f32], per_row: usize) {
    for row in values.chunks(per_row) {
        let _ = writeln!(out, "{}", reals(row));
    }
}
