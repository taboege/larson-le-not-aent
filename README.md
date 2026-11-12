# Larson-Le verification

One Python program recomputes six Ingleton-Main (IM) extensions of Larson-Le's
rank-four dual matroid, contracts one element, and determines all 1,024 ranks
of a ten-element restriction by linear programming. It then checks the matroid
rank axioms and a violation of planar Desargues. Together with the mathematical
extension and minor theorems described below, this certifies that the starting
matroid is not almost entropic.

All computation and example data are in [`im_projection.py`](im_projection.py).
The program needs no input files, saved certificates, SageMath, or Anaconda.
It follows a fixed sequence of six extensions; it does not repeat the search
that discovered that sequence.

## Install and run

Use Python 3.10 or newer with `pip` and `venv`. The two runtime dependencies are
NumPy and the HiGHS Python interface, `highspy`. Installation reads their versions
from [`pyproject.toml`](pyproject.toml) and installs them automatically. HiGHS is
pinned to 1.16.0, the solver version used for verification.

Download or copy this repository, then open a terminal in its directory.

On Linux or macOS:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install .
larson-le-verify
```

On Windows, using PowerShell:

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install .
.\.venv\Scripts\larson-le-verify.exe
```

The virtual environment keeps the dependencies local to this installation.
The first installation needs access to the Python package index; computation
after installation is offline.

With the environment's Python, these are equivalent ways to run the program:

```sh
python -m im_projection
python im_projection.py
```

The second form requires the script in the current directory. You can also copy
just that file to a computer where the declared NumPy and HiGHS versions are
already installed. Python's standard library alone is insufficient.

Assertions are part of the verification: run without `python -O` or
`PYTHONOPTIMIZE`. The command refuses to run if assertions are disabled.

## Output

The default output directory is `im_projection_output` in the current directory.
To choose another directory and also export the six LP models:

```sh
larson-le-verify --output results --write-lps
```

On Windows, use `.\.venv\Scripts\larson-le-verify.exe` for the command above.
Reusing an output directory overwrites files with the same names.

* `extensions.json` records the six extension triples, the verified premises,
  and exact rational bounds for the closed sets after each compression.
* `projection.csv` records all 1,024 lower and upper rank bounds after contracting
  element 13 and restricting to `(1, 2, 4, 5, 6, 8, 12, 14, 15, 16)`. Bit `i` of
  `local_mask` denotes position `i` in that list, counting from zero. The `points`
  column uses the original element labels.
* With `--write-lps`, `extension_11.lp` through `extension_16.lp` contain the
  extension models before compression. Their objectives are zero; the program
  changes objectives while finding rank bounds.

A successful run ends by reporting that all nine Desargues collinear triples
have rank two, while the three proposed axis points have rank three. It also
reports 1,775 LP solves, including feasibility checks. Failure of a verification
check terminates the program with a nonzero exit status.

## What is checked

Unknown ranks are real variables throughout. The program imposes polymatroid
inequalities, verifies each extension's premises, and only identifies sets or
fixes ranks when exact lower and upper bounds coincide. Intermediate fractional
possibilities are retained.

HiGHS finds floating-point solutions. The program reconstructs rational feasible
points and checks their inequalities exactly. It also reconstructs nonnegative
dual multipliers and uses exact rational arithmetic and column bounds to obtain
valid objective bounds. Equality is certified by matching exact bounds, without
rounding ranks to integers. The final restriction is checked to have integral
ranks and satisfy the matroid rank axioms.

The mathematical ingredients outside the computation are the almost-entropic
three-line extension theorem, closure under minors, and the planar Desargues
theorem. The LPs themselves impose Shannon inequalities, which are necessary
conditions for almost entropicity. Relevant papers are:

* Matt Larson and Tuong Le, [*Duals of algebraic matroids need not be algebraic*](https://arxiv.org/abs/2609.02834).
  The initial non-bases are those of their dual matroid.
* Geva Yashfe, [*On the recognition problem for limits of entropy functions*](https://arxiv.org/abs/2509.06302),
  Theorem 18 and Remark 19, for planar Desargues.

## Tested environment and resources

A fresh installation on Linux x86-64 with Python 3.13.9, NumPy 2.5.4, and
HiGHS 1.16.0 passed the complete verification and exported all six LP models.
On an Intel Core i5-1335U, that run took 73.6 seconds and peaked at 480 MiB of
resident memory. HiGHS uses one solver thread. Runtime varies with hardware.
These measurements include recomputing everything from the embedded initial
matroid and extension triples; no saved results were loaded.

## Packaging

This is a single-module distribution using the standard
[`pyproject.toml` format](https://packaging.python.org/en/latest/guides/writing-pyproject-toml/)
and setuptools. `pip install .` installs the module, its dependencies, and the
`larson-le-verify` command. No separate `setup.py` or dependency list is needed.

To build a wheel and source archive for distribution, run these commands in an
activated virtual environment:

```sh
python -m pip install build
python -m build
```

The files appear in `dist/`. Someone receiving the wheel can install it with
`python -m pip install path/to/larson_le_verification-0.1.0-py3-none-any.whl`;
pip will also install its dependencies.
