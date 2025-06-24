![image](_static/arsenal.jpg)

# The Arsenal Gear Population Synthesis Code

Arsenal Gear is a Stellar Population Synthesis (SPS) code. The main
driving goals of its development are:

1.  **Flexibility** of choices related to input physics
2.  **Ease of use** from both the command line and within a python
    script
3.  **Stellar Feedback** is our primary modeling goal
4.  **Binary Stars** are easily modeled within the code

Arsenal Gear is intended to work in modes where stars are both sampled
from an IMF (stochastic mode) and when it is meant to be representative
of a well-sampled population and therefore quantities are computed via
integrals over the IMF (integrating mode). It is also intended to work
both for calculations from within a piece of python code or iPython
notebook as well as being able to be called from the command line or
from within C++/Fortran code for compatability with large hydrodynamical
simulations.

This code has been a collaborative effort on the parts of Eric
Andersson, Claude Cournoyer-Cloutier, Ben Keller, Lachlan Lancaster,
Marta Reina-Campos, James Wadsley, and Zachary Wyatt.

## Contents

```{toctree}
:maxdepth: 2
quickstart.md
parameters.md
outputs.md
models.md
code-models.md
architecture.md
troubleshooting.md
tests.md
docs.md
api.md
```