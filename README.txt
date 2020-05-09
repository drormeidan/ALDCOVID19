The simulations accompanying the arXiv paper in https://arxiv.org/abs/2004.01453 were
written in Matlab. There are two files:
CM4.m		Contains the ODEs and solves them numerically (forward Euler method);
CovidModel4.m	Sets parameters, calls CM4.m and creates the various plots.

As a rule, variable names are identical in the paper and the code. Variables not
appearing in the paper are documneted in the beginning of CM4.m.

Parameters:
[p_NS,p_M,p_S,p_C,p_HR,p_HD,p_VR,p_VD] = [0.3,0.78,0.14,0.08,0.85,0.15,0.5,0.5] % transition probabilities
[P1,P2,P3,P4] = [4.41,11.04,3.31,5.52] %scale parameters of weibull functions.
[P1,P2,P3,P4] = [1.47,1.47,1.47,1.47]  %shape parameters of weibull functions.
