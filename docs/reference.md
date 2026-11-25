# Reference programs

`reference/` holds the Kikuchi and Kanamori inversion package from the ETAL
distribution of the Earthquake Research Institute, University of Tokyo
(<https://wwweic.eri.u-tokyo.ac.jp/ETAL/KIKUCHI/>). The manual,
`reference/manual.pdf`, says the package may be freely copied; work that uses it
acknowledges M. Kikuchi and H. Kanamori.

## Build and regenerate the goldens

```bash
make -C reference            # programs in reference/build/bin/
make -C reference golden     # reruns both examples into reference/golden/
```

The build needs gfortran. [reference/Makefile](../reference/Makefile) compiles
`subroutine/*.f` into `mylib.a` and `recipes/*.f` into `recipes.a`, then links
each program in `main/` against both.

[reference/regenerate-golden](../reference/regenerate-golden) runs each example
in a scratch directory and copies every file the programs wrote, plus their
console output as `stdout`, into `reference/golden/<case>/`. Two runs give
identical bytes. `--inputs DIR` takes the four waveform and structure files
from `DIR`, and `--out DIR` writes the cases under `DIR`; the parity test uses
both.

| Case | Input | Chain | Result in `stdout` |
| --- | --- | --- | --- |
| `denali` | `fort.1`, `fort.2` | `green`, `inversion`, `graphics`, `inversion2`, `graphics2`, `plotcm` | Mw 7.78, variance 0.2606 |
| `peru` | `fort.1.Peru`, `fort.2.Peru` | `green3`, `inversion3`, `graphics3`, `mom3`, `gr3` | Mw 8.25, maximum slip 5.33 m |

## Program chains

Each program reads and writes fixed `fort.N` names in the working directory.
Most open their parameter file `i_<program>` themselves; `graphics`,
`graphics2` and `inversion2` read it from stdin.

| Program | Reads | Writes |
| --- | --- | --- |
| `green` | `fort.1`, `fort.2`, `i_green` | `fort.3` (Green's functions, binary) |
| `inversion` | `fort.1`, `fort.3`, `i_inversion` | `fort.16` (subevents), `fort.4` |
| `graphics` | `fort.4`, stdin | `plot8` (PostScript) |
| `inversion2` | `fort.1`, `fort.2`, stdin | `fort.24`, `fort.26` |
| `graphics2` | `fort.24`, stdin | `plot28` (PostScript) |
| `plotcm` | `fort.1`, `fort.3`, `i_plotcm` | `plot_cm.ps` |
| `green3` | `fort.1`, `fort.2`, `i_green3` | `fort.33` (binary) |
| `inversion3` | `fort.1`, `fort.33`, `i_inversion3` | `fort.36`, `fort.34` |
| `graphics3` | `fort.34`, `i_graphics3` | `plot38.ps` |
| `mom3` | `fort.1`, `fort.2`, `i_mom3` | `fort.24`, `fort.37` (slip in space and time), `fort.46` |
| `gr3` | `fort.1`, `fort.24`, `i_gr3` | `plot3.ps` |
| `plotw` | `fort.1`, `i_plotw` | `plot_w.ps` |
| `gregra` | `fort.3`, `i_gregra` | `plot_g.ps` |

`conv.sac.farm` and `conv.sac.iris` turn SAC records into `fort.1`, using the
J-B tables `jb.ptime`, `jb.stime` and `jb.table`. `rotSH` rotates horizontal
records to SH. `lists` and `listg` print `fort.1` and `fort.3`. The scripts in
`reference/workflow/` (`mk_conv.farm.pl`, `rdseed2sac`, `sacmacro`) drove that
conversion with rdseed and SAC. [architecture.md](architecture.md#prep)
describes the port of `conv.sac.farm`.

## Input files

`fort.1` holds one block per station after a header line with the event: ID,
station name, then azimuth, back azimuth, distance, slowness, geometrical
factor and `ix0`; the `im ib ic` and `izp izz ip` flags; magnification and
normalisation; maximum amplitude, sample count and sample interval; then the
samples. `fort.2` holds the structure: an ID line, then `TQP`, `TQS` and the
source layers, the receiver layers, and the PP bounce layers. The manual gives
every field.

## Changes from the distribution

- `main/conv.sac.iris.f`: in `GEOM`, `IMPLICIT REAL*8` moved above the
  `PARAMETER` statement. gfortran rejects the original order.
- `recipes/` keeps only the four Numerical Recipes routines the programs call:
  `eigsrt`, `jacobi`, `lubksb` and `ludcmp`.
- The build uses `-fno-automatic -finit-local-zero`. The code reads locals
  before it sets them and expects them to keep their values between calls.
  Without the flags, `plotcm`, `mom3` and `gr3` give different output on each
  run.

## MATLAB scripts

`reference/matlab/` holds post-processing scripts by Cesar Jimenez, in
Spanish: `formato.m` converts `fort.46` to geographic coordinates,
`espacio_temp.m` plots the space-time slip from `fort.37`, `replicas.m` plots
aftershocks, and `interp_def.m` interpolates deformation at geodetic stations.
`fotogramas.m` (Bento Caldeira, 2004) draws rupture frames and
`pfalla2latex.m` turns a fault file into a LaTeX table.
