The simulations accompanying the arXiv paper in https://arxiv.org/abs/2004.01453 were
written in Matlab. There are two files:
CM4.m		Contains the ODEs and solves them numerically (forward Euler method);
CovidModel4.m	Sets parameters, calls CM4.m and creates the various plots.

As a rule, variable names are identical in the paper and the code. Variables not
appearing in the paper are documneted in the beginning of CM4.m.
