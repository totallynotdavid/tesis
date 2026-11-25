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

## Fetch

`kk fetch` and `kk prep` ([src/main.rs](../src/main.rs)) take an event file,
[docs/events.md](events.md). [src/event.rs](../src/event.rs) parses it into a
`Config`, which gives the `Event` and P `Window` that `prep` reads. `main`
does the HTTP and file I/O; [src/fetch/](../src/fetch/mod.rs) builds the
requests and reads the answers.

1. `fetch::station_url` asks fdsnws-station for the selected channels within
   the event's distance range, at `level=response`.
2. [stationxml.rs](../src/fetch/stationxml.rs) reads each channel's
   coordinates and its first poles-and-zeros stage. A stage written in Hz is
   converted to rad/s. A channel whose response cannot be used keeps the
   reason, and `prep` reports it.
3. `fetch::bulk` writes one fdsnws-dataselect line per channel, covering its
   P window with 60 s on each side.
4. [miniseed.rs](../src/fetch/miniseed.rs) decodes the answer with EarthScope's
   `libmseed` crate into one trace per channel.
5. `fetch::record` joins a channel to its trace. A channel with no trace, or
   with a trace in more than one segment, is skipped.

## Prep

[src/prep/](../src/prep/mod.rs) is `conv.sac.farm`: it turns one record into
one `fort.1` station block. `convert` takes an `Event`, the `Earth` for its
depth, a `Record` that carries its own poles and zeros, and a `Window`. It
returns a `Station` or the `Skip` that names why the record was dropped.

1. [earth.rs](../src/prep/earth.rs) reads `jb.ptime`, `jb.stime` and
   `jb.table`; `Tables::jb()` parses the copies built into the binary. `Tables::at(depth)` traces the rays for that depth.
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
