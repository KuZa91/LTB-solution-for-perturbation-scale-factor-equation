# LTB-solution-for-perturbation-scale-factor-equation
**Paolo Marcoccia<sup>1</sup>, Giovanni Montani<sup>2</sup>**

<sub>1. University of Stavanger, Institutt for Matematikk og Fysikk, Kjølv Egelands hus, 5.etg, E-blokk, 4021 Stavanger, Norway </sub>  
<sub>2. University "La Sapienza" Roma, Dipartimento di Fisica, Piazza A.Moro 5, 00185 Roma, Italy</sub>  

## Introduction ##

In the following, we will analyze a weakly deformed isotropic Universe toward
a spherically symmetric model (see [G. Montani et al.](https://inspirehep.net/literature/896693),[P. J. E. Peebles](https://inspirehep.net/literature/376248)), thought as the natural metric framework of the [ΛCDM](https://arxiv.org/abs/astro-ph/9805201) model. We
find that inhomogeneities will arise in the previous context as small perturbations of the order of
few point percent over the background _FRW_ universe. The obtained picture offers a useful scenario
to investigate the influence of the inhomogeneity spectrum (left free in the obtained solution), on
the propagation of photons or gravitational waves at low redshift values, and in line of principle
may be used to account for several present-stage cosmological problems, such as the [Hubble Tension](https://academic.oup.com/mnras/article-abstract/doi/10.1093/mnras/stz3094/5849454?redirectedFrom=fulltext).
We present in this directory two different softwares built in _Fortran 95_:

- The [First One](https://github.com/KuZa91/LTB-solution-for-perturbation-scale-factor-equation/blob/master/LTBMatterCostant.f90) will solve _eqt. (75)_ of the [paper](https://arxiv.org/abs/1808.01489v3) using an _explicit euler_ method, the output of the code will be a _.txt_ having all the informations of the simulated variables needed to generate _figs. (1),(2),(3)_ ; 

- The [Second One](https://github.com/KuZa91/LTB-solution-for-perturbation-scale-factor-equation/blob/master/MinLambdaEstimate.f90) will instead find the zeros of _eqt. (75)_ by using a _divide et impera_ method, the output would be a _.txt_ having the estimated value of &Omega; &Lambda; on each iteration step, and may be used to generate _fig. (4)_;

## Analysis Details ##

Details of the analaysis can be found in our [preprint paper](https://arxiv.org/abs/1808.01489v3).

## Additional Material ##

- In order to reproduce the results of our paper, a simple installation of [fortran 95](https://gcc.gnu.org/wiki/GFortran) is needed;

- The additional files [.plt](https://github.com/KuZa91/LTB-solution-for-perturbation-scale-factor-equation/blob/master/LTBMatterCostanta.plt) will be used by the main script for generating a first plot of the simulated variables, and require [gnuplot](http://www.gnuplot.info/) in order to be run;

- The directory, contains some of the [figures](https://github.com/KuZa91/LTB-solution-for-perturbation-scale-factor-equation/blob/master/EvolutionOfA.png) used for the [paper](https://arxiv.org/abs/1808.01489v3), these figures were generated using the [Anaconda](https://www.anaconda.com/) distribution of python, as well as the [seaborn](https://seaborn.pydata.org/) data visualization library;

- The fits shown in _figs. (5),(6)_, where done by using the [OriginLab](https://www.originlab.com/) software, additional information on the fits used are reported in the [paper](https://arxiv.org/abs/1808.01489v3).


## Additional information about the execution of the analysis

- In order to avoid errors, tt is highly reccomended to delete the first comment rows of the _.txt_ file before loading that on third part softwares as [Pandas](https://pandas.pydata.org/) or [OriginLab](https://www.originlab.com/);
                                     
- If you wish to reproduce the same plots of the [paper](https://arxiv.org/abs/1808.01489v3), as done for [EvolutionOfA.png](https://github.com/KuZa91/LTB-solution-for-perturbation-scale-factor-equation/blob/master/EvolutionOfA.png), a guide on how to generate fancy plots may be found in the [Jake VanderPlas D.S Handbook](https://jakevdp.github.io/PythonDataScienceHandbook/04.08-multiple-subplots.html).
