# Solving a Search-Theoretic Model of Money w/ Persistent Heterogeneity

By: Christian Bustamante <br>
<a href="https://cbustamante.co"><https://cbustamante.co/></a><br>


## About

This set of codes solves a search-theoretic model of money similar in spirit to Lagos and Wright (2005), but where
preferences are not quasi-linear. Dropping the quasi-linearity assumption generates a distribution of
money that becomes persistent across periods. In this model, agents accumulate only one asset (money).

For a description of the model and the numerical method used to solve it, see the appendix for
"*The Long-Run Redistributive Effects of Monetary Policy*" in [my research webpage](https://cbustamante.co/research).
The file `apx_1asset.pdf` in this repo, provides the relevant pages of said appendix.

## Contents

This package is organized as follows:

- The folder `src` contains the main codes that solve the model. The main program is in `Main_Only_Money.f90`.
- The folder `lib` provides a set of different numerical subroutines needed by the codes in `src`. 
- The folder `plot` contains Matlab codes to plot/analyze some of the model results. These codes reproduce the figures in Appendix D.
- In the root folder, there is the `Makefile` used to build the solution program. It is set up to run with `ifort` but it should work (if adapted) with other compilers. It also uses the `hdf5`, `mkl`, and `openmp` libraries.

## License and Citation

MIT license. <br>
Citation: <br>
Bustamante, C. (2023). "The Long-Run Redistributive Effects of Monetary Policy," *Journal of Monetary Economics*, Forthcoming.
