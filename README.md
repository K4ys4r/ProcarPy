# ProcarPy Module
## What is ProcarPy
**ProcarPy** is a [**Python**](https://www.python.org/) module which allows to parse and plot the electronic band structure diagram from [**VASP**](https://www.vasp.at/) PROCAR file.<br>
[**Numpy**](https://www.numpy.org/) and [**Matplotlib**](https://matplotlib.org/) packages are required.

## What ProcarPy can do ?
There is different [**PROCAR output format**](https://cms.mpi.univie.ac.at/wiki/index.php/PROCAR). However ProcarPy can parse and process all PROCAR files from tags below:

- LORBIT = 10 : PROCAR not decomposed (_s_, _p_, _d_).
- LORBIT = 11 : PROCAR lm decomposed (_s_, _p<sub>y</sub>_, _p<sub>z</sub>_, _p<sub>x</sub>_, _d<sub>xy</sub>_, _d<sub>yz</sub>_, ...).
- ISPIN = 2 : Spin Polarized Calculation (_s<sub>up</sub>_, _s<sub>down</sub>_, _p<sub>up</sub>_, _p<sub>down</sub>_, ...).
- LNONCOLLINEAR=.TRUE. : Spin-orbit coupling.<br>

Then, with ProcarPy, the total and projected electronic band structure can be plotted.

## Usage of ProcarPy
Some tutorials have been performed to show how ProcarPy can be handled.<br>
- [**Si Diamond Structure**](https://github.com/K4ys4r/ProcarPy/blob/master/test/Si_diamond/Si_Diamond_Tutorial.ipynb)
- [**Hematite (&alpha;-Fe<sub>2</sub>O<sub>3</sub>) Spin-Polarized calculation**](https://github.com/K4ys4r/ProcarPy/blob/master/test/Fe2O3/Fe2O3_Tutorial.ipynb)
- [**Cobalt Spin-Orbit coupling**](https://github.com/K4ys4r/ProcarPy/blob/master/test/Co/Co_Tutorial.ipynb)
