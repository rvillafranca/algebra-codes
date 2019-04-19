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
high level capabilities to automate the generation of the code basis/generator matrix. Furthermore the Python/Sage routines can split the code for parallel
processing.

The usage of these routines can be understood by looking on the comments in the
source files and the example programs. For those who intend to use the C
routines directly, it should be noted that there are two different routines, one
for codes over prime fields and other for codes over non-prime fields.

**TODO**: Push the Python/Sage stuff. (Coming really soon.)
