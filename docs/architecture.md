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

## Prep

[src/prep/](../src/prep/mod.rs) is `conv.sac.farm`: it turns one record into
one `fort.1` station block. `convert` takes an `Event`, the `Earth` for its
depth, a `Record` that carries its own poles and zeros, and a `Window`. It
returns a `Station` or the `Skip` that names why the record was dropped.

1. [earth.rs](../src/prep/earth.rs) reads `jb.ptime`, `jb.stime` and
   `jb.table` once. `Tables::at(depth)` traces the rays for that depth.
   `Earth::geometry` gives azimuths, distance, slowness and spreading.
   `Earth::travel` gives the P or S onset.
2. `convert` cuts `duration` seconds from `lead` seconds before the onset and
   removes the mean.
3. It divides out the response, band-passes with an eighth-order zero-phase
   filter, and converts to the event's quantity in the frequency domain. The
   FFT and the complex arithmetic are in [fft.rs](../src/prep/fft.rs).
4. It decimates to the window's interval and normalises by the peak.

`convert` differs from `conv.sac.farm` in four places:

- The station's `power` is the number of times the output quantity
  differentiates displacement, which is how `green` reads it. The reference
  writes 1 for velocity and 0 otherwise, so `green` computes displacement
  synthetics for acceleration records.
- Samples keep full precision. The reference prints three decimals.
- An output interval shorter than the record's is refused as
  `Skip::Upsampling`. The reference interpolates.
- Counts without deconvolution (`id = 3`) are not offered.

## Parity

[tests/parity.rs](../tests/parity.rs) parses the four example inputs, writes
them back with the Rust writers, and runs
[reference/regenerate-golden](../reference/regenerate-golden) on the result.
Every output file must equal its golden byte for byte. A change of 0.01 km/s in
one layer velocity changes all nine Denali outputs, so a writer that loses
precision fails the test.

[tests/prep.rs](../tests/prep.rs) builds a copy of `conv.sac.farm` whose
`fort.22` formats print nine significant digits, enough to identify an `f32`.
It writes six synthetic SAC records and a station inside 30 degrees, runs
the copy and `convert` on them for three output quantities and depths, and
requires every value to be equal. A one-ulp change in pi fails it.
