# Figures from Lizier et al. 2023

This folder contains scripts / parameters files to recreate the results of:

J.T. Lizier, F. Bauer, F.M. Atay, and J. Jost,
_"Analytic relationship of relative synchronizability to network structure and motifs"_,
Under submission,
2023

## Figure 1 - convergence of numerical results from simulations

The parameter settings for recreating figure 1 are contained in `figure1Parameters.m`.

The script to run the experiments is `figure1ConvergenceOfSimulations.m`.
Run that script in Matlab from this folder:
```matlab
figure1ConvergenceOfSimulations
```

An `.eps` and `.mat` file is saved for the figure.

*Crucially* - note that the script adjusts the number of repeat runs 
downwards from 2000 in the parameters file (matching the paper) to 10 -- this
is done to ensure that the results run quickly for you. 
Also, the script removes the largest number of time series samples (L=1,000,000)
for the same reason.
These shorter runs take about 8 minutes 
on my machine.

Obviously you can remove those restrictions -- running for 2000 repeats and with more time series samples
will take about 100 hours (the longest number of time series samples really 
blows this out)
so you would be advised to run for more repeats in parallel fashion as per the
[cluster instructions](/cluster).

In any case, the shorter runs give a good idea of the main trends already,
though there are some fluctuations / std error differences compared to the longer runs.

## Figure 4 - small world sweeps

The parameter settings for recreating figure 4a are contained in `figure4aParameters.m`,
and these parameters are then adjusted by the experimental script for the other
sub-figures.

The script to run the experiments is `figure4SyncThroughSmallWorld.m`.
Run that script in Matlab from this folder:
```matlab
figure4SyncThroughSmallWorld
```

Notice that this script recreates each sub-figure in order, and modifies
the parameters object loaded in from `figure4aParameters.m` to have the correct
parameters for each sub-figure.
An `.eps` and `.mat` file is saved for each sub-figure.

*Crucially* - note that the script adjusts the number of repeat runs 
downwards from 1000 in the parameters file (matching the paper) to 10 -- this
is done to ensure that the results run quickly for you. 
The run with 10 repeats for each parameter combination takes about 12 minutes 
on my machine.

Obviously you can remove those restrictions -- running for 1000 repeats will take about 20 hours,
so you would be advised to run for more repeats in parallel fashion as per the
[cluster instructions](/cluster).

In any case, running 10 repeats only gives a good idea of the main trends already,
though there are some fluctuations / std error differences compared to the longer runs.