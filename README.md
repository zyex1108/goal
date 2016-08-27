# Goal
> Adaptive goal-oriented solid mechanics

Goal is a mini-application written in C++11
that implements goal-oriented error estimation
and adaptation for solid mechanics applications.
The main emphasis of Goal is to provide reliable
and accurate error estimation for physically
meaningful functional quantities, and the
ability to reduce these errors through mesh
adaptation.

## Installation / Getting Started

### Dependencies
To get started, you need
* [Cmake][0] greater than version 3.0.0,
* A valid C++11 compiler with [MPI][1] wrappers
* A valid [PUMI][2] installation,
* A valid [Trilinos][3] installation.

For detailed dependency build instructions,
see the [build instructions][4] wiki page.

### Configuration
An example configuration script for Goal
is located in [here][5] in the repository.

#### GOAL_FAD_SIZE
Default: `16`

Goal makes use of a process called forward
[automatic differentiation][6]. This process
stores an array of derivative values for
each variable computed. This configuration
variable specifies the maximum size of the
derivative array, which should correspond
to the maximum number of elemental degrees
of freedom that will be used with a given
goal build. For example, a tetrahedral mesh
with linearly interpolated displacement and
pressure degrees of freedom will have
16 elemental degrees of freedom. Problems
with greater than GOAL_FAD_SIZE elemental
degrees of freedom will fail.

#### GOAL_DISABLE_CHECKS
Default: `OFF`

Goal performs basic sanity checks, such as
checking that user-specified boundary conditions
are applied to geometric boundaries that
actually exist. For extreme production runs,
where performance is absolutely critical, these
sanity checks can be turned off.

#### GOAL_TESTING
Default: `OFF`

A suite of simpler regression [tests][7]
is provided with this repository for
developers and to check build validity.
This configuration option enables the
ability to run these tests under CTest.

[0]: https://cmake.org
[1]: https://www.mpich.org
[2]: https://github.com/scorec/core
[3]: https://github.com/trilinos/Trilinos
[4]: https://github.com/bgranzow/goal/wiki/Build-Instructions
[5]: https://github.com/bgranzow/goal/blob/master/aux/do-config.sh
[6]: https://en.wikipedia.org/wiki/Automatic_differentiation
[7]: https://github.com/bgranzow/goal/tree/master/tests
