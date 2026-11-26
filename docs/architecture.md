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
   FFT and the complex arithmetic, shared with `green`, are in
   [src/fft.rs](../src/fft.rs).
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

The reference leaves bad channels to the analyst. `outliers` drops them once
every record is converted: a station whose peak over its spreading factor is
more than `OUTLIER` (5) times the median of its component's stations is
skipped as `Skip::Amplitude`. Yauca's good stations span 0.7 to 1.5 times the
median; IU.RCBR, whose borehole channel drifts with no zero crossing, is 22.

## Green

`kk green` reads `fort.1` and the event file's structure and writes `fort.3`.
[src/green/](../src/green/mod.rs) is `green` and `green3`. `base` gives the
six base moment-tensor elements per station and source point, for `fort.3`.
`fixed` gives one mechanism, for `fort.33`. Each synthetic is built in the
frequency domain and transformed back with `cfft`:

1. [response.rs](../src/green/response.rs) gives Q (`qfm`) and the instrument,
   from poles and zeros (`instg`) or the WWSSN constants (`instw`).
2. [radiation.rs](../src/green/radiation.rs) is `radp`: the P, SV and SH
   radiation of a mechanism, up- and down-going.
3. [propagator.rs](../src/green/propagator.rs) applies Haskell's layer
   matrices: `refl` and `reflsh` give the source-side up- and down-going
   waves for the depth phases, `cnvr` and `cnvrsh` convert the incident wave
   to ground motion under the receiver, and `refl` gives the PP bounce
   point.
4. `Model::bodyw` combines them as `sub.bodyw.f` does.

The port mirrors gfortran's `-O` arithmetic to reach bit parity: a real
operand in a complex expression is widened to complex first, products and
quotients are full complex operations with Smith's division, and `x**3` is
`(x*x)*x`.

`bodyw` differs from the reference in one place. The reference computes the
P delay `tp` for P and PP only, so an SV synthetic uses the `tp` its previous
call saved (`-fno-automatic`). `bodyw` gives SV its own.

[src/formats/greens.rs](../src/formats/greens.rs) and
[fixed.rs](../src/formats/fixed.rs) are `fort.3` and `fort.33`, gfortran
unformatted files; [unformatted.rs](../src/formats/unformatted.rs) reads and
writes their records.

## Inversion

`kk invert` reads `fort.1`, `fort.3` and the event file's `[inversion]` table
and writes `fort.16` and `fort.4`. `Config::parameters` in
[src/event.rs](../src/event.rs) turns the table into the `i_inversion`
[Parameters](../src/formats/inversion.rs). It expands the station-code weights
to one per station `green` kept, in `fort.1` order.

[src/inversion/](../src/inversion/mod.rs) is `inversion`. `invert` runs these
stages:

1. `stime` builds the source time function. `records` convolves it with each
   station's base elements and gives the station's delay at every grid point
   and depth.
2. `inverse` solves the normal equations at each depth with `ludcmp` and
   `lubksb`.
3. `search` finds each subevent: the grid point, onset and weights that most
   reduce the residual. With more than one subevent, two more rounds refit
   each one against the records with the others removed.
4. `solution` sums the subevents and gives the standard deviation and the
   moment magnitude.

[linalg.rs](../src/inversion/linalg.rs) holds the Numerical Recipes routines
and the cubic and quadratic solvers. [mechanism.rs](../src/inversion/mechanism.rs)
holds the determinants, the double-couple constraint `zmom`, and `dcp`, which
splits a tensor into a double couple and its axes.
[src/formats/solution.rs](../src/formats/solution.rs) writes the `Solution`
as `fort.16` and the summary line, and [plot.rs](../src/formats/plot.rs) as
`fort.4`. [fortran.rs](../src/formats/fortran.rs)
formats reals and integers as gfortran's `F` and `I` descriptors do.

Fortran `REAL` is `f32`, so the port computes in `f32`. It keeps the
reference's constants, its order of operations, and comparisons that send NaN
down the same branch. The module allows the clippy lints those choices trip.

`invert` differs from the reference where the reference reads past an array
or runs on bad input:

- A phase outside 1 to 4, a first round with no fitting grid point, singular
  normal equations, and a tensor whose eigenvalues do not converge are
  errors. The reference runs on with stale values or pauses.
- An onset that puts a correlation lag outside the record reads 0.
- A source-time flank longer than the record is cut at the record.
- A reference point `k0` outside the line of depths gives zero weights for
  the standard deviation.
- There are no array limits. The reference holds 1024 samples, 11 grid points
  and 11 depths, 50 subevents and 40 stations.
- The reference calls `dcp` at every improved grid point and discards the
  result. `invert` calls it once per subevent.

## Plot

`kk plot` reads `fort.4` into a [Plot](../src/formats/plot.rs) and draws it as
two SVG files. It replaces `graphics`, which writes PostScript.
[src/plot/mod.rs](../src/plot/mod.rs) lays out the figures and
[svg.rs](../src/plot/svg.rs) writes the elements; there is no drawing crate.

[sphere.rs](../src/plot/sphere.rs) holds the focal sphere in the reference's
frame: x north, y east, z down, on a lower-hemisphere equal-area net. A ball
fills where the tensor's P amplitude `dᵀMd` is positive, row by row, and draws
the nodal planes from strike, dip and slip with Aki and Richards' normal and
slip vectors. Stations sit at their azimuth and the takeoff angle
`asin(p v)`, where `fort.4` stores the slowness times the source velocity. [tests/plot.rs](../tests/plot.rs) checks on the
Denali `fort.4` that each nodal plane lies where its tensor radiates nothing
and that the T axis is compressional.

## Sweep

`kk sweep` ([src/main.rs](../src/main.rs)) reads `fort.1` once and runs the
Green's functions, the inversion and the figures for each structure file.
`Files::under` in [src/event.rs](../src/event.rs) moves the outputs into a
directory per structure. [src/compare.rs](../src/compare.rs) gives Kagan's
angle between the total and the event file's `[catalogue]` double couple,
from the T, P and B axes of each.

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

[tests/green.rs](../tests/green.rs) requires `base` on the Denali example
and `fixed` on the Peru example to write the golden `fort.3` and `fort.33`
byte for byte. The goldens hold only P vertical and SH records on GDSN
instruments without poles, so the test also builds `green` with an SV `tp`
patched in. It runs that build on Denali stations rewritten to SV, PP,
horizontal components, WWSSN and poles and zeros, over four depths, under two
structures. Every sample must match.

[tests/inversion.rs](../tests/inversion.rs) requires `invert` on the Denali
example to write the golden `fort.16`, `fort.4` and summary line byte for
byte. The golden runs one basis, one time function and one depth, so the test
also builds the reference `inversion` and runs both on Denali over six cases.
The cases cover every basis and time function, points and lines, one and three
depths, one to thirteen subevents, and stations left out by weight 0.
