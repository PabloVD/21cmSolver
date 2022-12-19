# 21cmSolver

[![arXiv](https://img.shields.io/badge/arXiv-1912.09488-B31B1B.svg)](http://arxiv.org/abs/1912.09488)
[![arXiv](https://img.shields.io/badge/arXiv-2004.00013-B31B1B.svg)](http://arxiv.org/abs/2004.00013)

Python codes to compute the 21 cm cosmological global signal. It allows to evolve the signal for a grid of astrophysical parameters.

It allows to take into account the contribution from non-linear cosmological overdensities through an analytical framework, as developed in [Phys. Rev. D 101, 083502](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.101.083502) ([arXiv:1912.09488](https://arxiv.org/abs/1912.09488)). Contribution from the Ly-alpha radiative transfer is included as explained in
[JCAP, 2006(06):026](https://iopscience.iop.org/article/10.1088/1475-7516/2020/06/026) ([arXiv:2004.00013](https://arxiv.org/abs/2004.00013)). Please refer to these articles for more info.


## Code Description

Here is a description of the different programs:

- `21cmEvolver.py`: runs a simulation, solving the differential equations and computing the relevant quantities, such as the 21 cm brightness temperature.

- `Runner.py`: run the solver code `21cmEvolver.py` on a grid of astrophysical parameters, for a list of values of the initial overdensity Delta.

- `Plotter.py`: code to plot features such as the ionization fraction, temperatures, global signal and Lyman alpha flux.

In the `Source` folder:

- `constants.py`: some relevant physical and cosmological constants.

- `astro_params.py`: the relevant astrophysical parameters are included here. You may want to change these values (or run `Runner.py` to solve in a grid).

- `functions.py`: utility functions for solving the differential equations.

- `lyalpha_rt.py`: functions related to the Lyman-alpha radiative transfer.


## Requisites

This code makes use of several Python libraries:

* `numpy`
* `matplotlib`
* `scipy`
* [mpmath](https://mpmath.org/doc/current/index.html)
* [colossus](https://bdiemer.bitbucket.io/colossus/)


## Citation

If you use the code, please link this repository, and cite [arXiv:1912.09488](https://arxiv.org/abs/1912.09488) and/or [arXiv:2004.00013](https://arxiv.org/abs/2004.00013).


## Contact

Feel free to contact me at <pablo.villanueva.domingo@gmail.com> for comments, questions and suggestions.
