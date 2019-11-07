# 21cmSolver

Python code to compute de 21 cm global signal. It allow to evolve the signal for a grid of astrophysical parameters, taking into account over densities.
For more info, see arxivXXXX

Description of the different programs:

- 21cmEvolver.py: runs a simulation, solving the differential equations and computing the relevant quantities, such as the 21 cm brightness temperature.

- Runner.py: run the solver code 21cmEvolver.py on a grid of astrophysical parameters, for a list of values of the initial overdensity Delta.

- Plotter.py: code to plot the ionization fraction, temperatures, global signal and Lyman alpha flux.

In the Source folder: the needed functions and constants are included.

- constants.py: some relevant physical constants

- astro_params.py: the relevant astrophysical parameters are included here. You may want to change these values (it can be done in Runner.py)

- functions.py: functions for solving the differential equations

- lyalpha_rt.py: functions related with the Lyman-alpha radiative transfer

This code makes use of several Python libraries, such as numpy, scipy and colossus.
