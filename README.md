# Conditional_Constraint_Fitting
This example code demonstrates how to apply the conditional constraint method to fit a model prediction to a measurement. Intended as supplemental material to <arxiv>.

Model prediction with covariance matrix spanning measurement and parameter space is intended as input. This can be built from a set of universes, such as those generated in NUISANCE, or can be supplied directly. The (separate) applications of the regularization matrix and the quantile mapping mitigation strategy are demonstrated in example fits to cross sections from 10.1103/PhysRevLett.133.041801.

Build the example .cxx script with "make" and then run it. Included are fits with and without the regularization matrix applied, as well as with and without quantile mapping performed. Additionally, a basic framework is provided to fit a subset of all available model parameters. Note the need to reconfigure the class after applying the regularization matrix as this is introduced at an early stage of setup when MC files are read in.
