The simulations accompanying the arXiv paper in https://arxiv.org/abs/2004.01453 were
written in Matlab. 

The codes was seperated to three folders:
Model: Contain producing codes of figures relate to the model without Data Analysis.
	covidmodel5.m:   
		Producing Main Figure, PWQ and panel (n) and (o) in Fig.3.
	
	covidmodel5_defectors.m:
		Producing Defectors figure.

	covidmodel5_spi.m:
		Producing Selective Isolation figure.

	covidmodel5_tests.m:
		Producing Test Figure.

	BA_generator_fast.m:
		Generate ER or SF network with households cliques.

	temporalmatrix.m:
		Generate the temporal links of the network.

	closeS.m/closeS_spi.m:
		Isolated houeholds and filter the outdoor temporal links relative to quarantine strategies. 

DataAnalysis:
	estBeta.m:
		Estimate the beta parameter and produciong DataAnalysis figure. 

Grid:
	Grid_TinTout_AlphaBeta.m:
		Calculate grid of alpha and beta where given T_in and T_out values (with option to plots).
	closeS_grid.m:
		closeS version to grid.

Notes:
As a rule, variable names are identical in the paper and the code. Variables not
appearing in the paper are documneted in comment in their first instance.

The codes need the machine learning and statistics toolbox of Matlab.
Some of them used the parallel computing toolbox.

