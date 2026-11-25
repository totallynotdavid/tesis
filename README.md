# tesis

Source inversion of the Yauca (Arequipa, Peru) earthquake with the Kikuchi and
Kanamori teleseismic body-wave method. The repository holds the thesis text and
the reference Fortran programs, with golden outputs from two example events.

```bash
make -C reference golden
```

That builds the 18 reference programs with gfortran and reruns the Denali 2002
and Peru 2001 examples into `reference/golden/`.

- `reference/`: the Kikuchi-Kanamori Fortran 77 programs, example inputs,
  golden outputs, and the advisor's MATLAB post-processing scripts.
- `paper/`: the thesis proposal and lab reports, in LaTeX.

The reference build, the program chains and the file formats are described in
[docs/reference.md](docs/reference.md).
