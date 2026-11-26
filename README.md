# tesis

Source inversion of the Yauca (Arequipa, Peru) earthquake with the Kikuchi and
Kanamori teleseismic body-wave method. The repository holds the thesis text,
the reference Fortran programs with golden outputs from two example events, and
`kk`, a Rust port of the tooling checked against those outputs.

```bash
mise install
mise run check
```

That pins Rust, lints the crate, and runs the tests. The parity test needs
gfortran: it builds the reference programs and runs them on files the Rust
writers produced.

```bash
cargo run -- fetch events/yauca-2013/event.toml   # FDSN stations and records
cargo run -- prep events/yauca-2013/event.toml    # fort.1
cargo run -- green events/yauca-2013/event.toml   # fort.3
cargo run -- invert events/yauca-2013/event.toml  # fort.16 and fort.4
cargo run -- plot events/yauca-2013/event.toml    # solution.svg and waveforms.svg
```

[docs/events.md](docs/events.md) describes the event file.

- `src/`: the `kk` crate. It fetches FDSN records, converts them to `fort.1`,
  computes the Green's functions into `fort.3`, and inverts them for the
  subevents into `fort.16` and `fort.4`, and draws the solution as SVG.
- `events/`: event files, one directory per event.
- `reference/`: the Kikuchi-Kanamori Fortran 77 programs, example inputs,
  golden outputs, and the advisor's MATLAB post-processing scripts.
- `paper/`: the thesis proposal and lab reports, in LaTeX.

The code map is in [docs/architecture.md](docs/architecture.md). The reference
build, the program chains and the file formats are in
[docs/reference.md](docs/reference.md).
