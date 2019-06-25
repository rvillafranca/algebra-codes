# algebra-codes

Some routines are provided to compute weight distributions of linear codes over
finite fields.

The determination of weights is done by exhaustive computation of
the codewords. In order to improve efficiency this exhaustive computation is
done by C routines which generate the vector-space elements in such an order
that each new element is obtained from the previous one by a simple vector
addition (that is, some integer additions).

A Python wrapper is provided for these routines to be used by a higher level
procedure in [SageMath](https://www.sagemath.org/). This offer the possibility
of using Sage's resources for working with algebraic structures and Python's
high level capabilities to automate the generation of the code basis/generator
matrix.

The usage of these routines can be understood by looking on the comments in the
source files and the example programs. For those who intend to use the C
routines directly, it should be noted that there are two different routines, one
for codes over prime fields and other for codes over non-prime fields.

There is a minimal `Makefile` provided which uses `gcc`, which sholud be present
on any Linux distribution as well as on the SageMath environment on Windows (v.
remark below).

Running
```bash
make shared
```
will build in the `sage/` subfolder the shared library, `exhaustive.so`, to be
used by the SageMath routines.

Running
```bash
make examples
```
will build an executable binary that can be run in a terminal.

> **Remark**: There is a native Windows binary installer since SageMath 8.0
> (https://wiki.sagemath.org/SageWindows,
> http://www.sagemath.org/download-windows.html). It is based on Cygwin and
> comes with `gcc`. Among other ways, a UNIX-like environment can be accessed
> via Python's `os` package. To build the shared library, for example, one can
> do, from the folder where the `Makefile` is:
>
> ```python
> import os
> os.system("make shared")
> ```
