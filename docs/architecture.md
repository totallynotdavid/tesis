# Architecture

`kk` is one Rust package. Its modules replace the Fortran tooling in
`reference/` one stage at a time, and each stage is checked against the golden
outputs that [reference.md](reference.md) describes.

## Formats

[src/formats/](../src/formats/mod.rs) reads and writes the two inputs of every
program chain.

- [records.rs](../src/formats/records.rs) is Fortran list-directed input. A
  `Records` hands out one `Read` per Fortran `READ` statement. The read pulls
  values across records and drops the rest of the record it stops in. Repeat
  counts (`3*0.0`) and `/` are refused with their line number; no example file
  uses them.
- [waveforms.rs](../src/formats/waveforms.rs) is `fort.1`: a title, then per
  station its geometry, a WWSSN or GDSN instrument, and the samples.
- [structure.rs](../src/formats/structure.rs) is `fort.2`: a title, the
  attenuation times, and the source, receiver and bounce-point layers. The
  programs read only the first model in the file, and so does `parse`.

A writer emits each real as Rust's shortest round-trip decimal. gfortran reads
it back to the same bits.

## Parity

[tests/parity.rs](../tests/parity.rs) parses the four example inputs, writes
them back with the Rust writers, and runs
[reference/regenerate-golden](../reference/regenerate-golden) on the result.
Every output file must equal its golden byte for byte. A change of 0.01 km/s in
one layer velocity changes all nine Denali outputs, so a writer that loses
precision fails the test.
