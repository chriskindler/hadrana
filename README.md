# hadrana

Analysis pipeline for extracting nucleon observables such as masses, energies at
finite momentum, and the axial and pseudoscalar form factors from correlation function measurements
on CLS gauge ensembles. Requires [`chigrad`](https://github.com/chriskindler/chigrad),
which provides an implementation of [`iminuit`](https://github.com/scikit-hep/iminuit), the Python interface to [MINUIT2](https://root.cern/doc/master/Minuit2Page.html), which in turn is a C++ reimplementation of the original [MINUIT](https://seal.web.cern.ch/seal/snapshot/work-packages/mathlibs/minuit/)
minimiser of the theory divison at CERN.
