# 21cmSolver

Python code to compute the 21 cm global signal. It allows to evolve the signal for a grid of astrophysical parameters, taking into account the overdensities.
For more info, see [arXiv:1912.09488](https://arxiv.org/abs/1912.09488), [Phys. Rev. D 101, 083502](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.101.083502)

Description of the different programs:

- 21cmEvolver.py: runs a simulation, solving the differential equations and computing the relevant quantities, such as the 21 cm brightness temperature.

- Runner.py: run the solver code 21cmEvolver.py on a grid of astrophysical parameters, for a list of values of the initial overdensity Delta.

- Plotter.py: code to plot the ionization fraction, temperatures, global signal and Lyman alpha flux.

In the Source folder:

- constants.py: some relevant physical constants.

- astro_params.py: the relevant astrophysical parameters are included here. You may want to change these values (it can be done in Runner.py).

- functions.py: functions for solving the differential equations.

- lyalpha_rt.py: functions related with the Lyman-alpha radiative transfer.

This code makes use of several Python libraries, such as numpy, scipy and [colossus](https://bdiemer.bitbucket.io/colossus/).
