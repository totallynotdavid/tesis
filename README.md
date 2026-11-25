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

- `src/`: the `kk` crate. It reads and writes the `fort.1` waveform and
  `fort.2` structure files.
- `reference/`: the Kikuchi-Kanamori Fortran 77 programs, example inputs,
  golden outputs, and the advisor's MATLAB post-processing scripts.
- `paper/`: the thesis proposal and lab reports, in LaTeX.

The code map is in [docs/architecture.md](docs/architecture.md). The reference
build, the program chains and the file formats are in
[docs/reference.md](docs/reference.md).
