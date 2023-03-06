# Linear Sync Toolkit (linsync)

Copyright (C) 2012- [Joseph T. Lizier](http://lizier.me/joseph/)

The `linsync` toolkit provides Matlab tools for analysing **synchronization** in networks of linearly coupled nodes.
Specifically, it provides implementations of the maths to measure and explore
the expected mean square deviation from synchronization <\sigma^2> as a function of network coupling matrix C, with
the mathematical details provided in the following paper:

J.T. Lizier, F. Bauer, F.M. Atay, and J. Jost,
"Analytic relationship of relative synchronizability to network structure and motifs",
Under submission,
2023

The tools compute \sigma^2 in both:
* Continuous-time Ornstein-Uhlenbeck dynamics: dX(t) = -X(t)*(I-C) + dw(t), where
  w(t) is a multivariate Wiener process (uncorrelated) with covariance matrix I,
  and C is the update matrix, and
* Discrete-time auto-regressive processes: X(t+1) = X(t)*C + R(t), where
  where R(t) is uncorrelated mean-zero unit-variance Gaussian noise
 and C is the update matrix.

In both cases the NxN weighted connectivity matrix C above is specified in row-vector form,
the mathematics handles the case where the fully synchronized state vector
\psi_0=[1,1,...,1] is an eigenvalue of C with \lambda_0 = 1,
and <\sigma^2> = \lim_{t -> \infty} < 1/N \sum_i (x_i(t) - \bar{x(t)} )^2 >

Please *cite* your use of the toolkit via the above paper.

# Use cases

In this section we briefly outline the primary use cases here are:
1. Generating sample network sructures _C_ to investigate, and
1. Computing the deviation from sync \sigma^2 for a given network structure _C_.
Then building on those we discuss more involved use cases for:
1. Batch experiments involving parameter sweeps and repeat runs, or
1. Running experiments on a cluster.
Finally, with the results generated you can:
1. Plot the results from these runs.

## Generating network structure

In order to compute <\sigma^2> for a network, we need a directed weighted connectivity matrix _C_ for it.
Several scripts are distributed for the user to generate _C_ matrices for standard structures, including:
* Erdos-Renyi random graphs -- `generateNewRandomMatrix.m`
* Fixed in-degree random graphs -- `generateNewRandomMatrixFixedD.m`
* Watts-Strogatz ring networks with fixed in-degree -- `generateNewRandomRingMatrix.m`

These scripts generate unweighted adjacency matrices, then the user should weight the edges.
A sample generating a weighted, directed Watts-Strogatz ring network of _N_ nodes, with in-degree _d_,
and rewiring probability _p_, without self-connections allowed from re-wiring, and ensuring 
the network remains connected:

```matlab
% First generate the unweighted adjacency matrix without self-edges:
A = generateNewRandomRingMatrix(N, d, p, false, undirected, false, true);
% Then compute the in-degrees (should be all d in this example):
D = diag(sum(A));
% and weight the edges as c/d each and insert self-edges with weight b - c:
% (note that using b=1 ensures \psi_0 is an eigenvector of C with \lambda_0 = 1)
C = (b - c) .* I + c .* A * inv(D);
```

## Use cases of analysing the sync for a given network _C_

Computing the expected deviation <\sigma^2> from the synchronized state for
a given network connectivty matrix _C_ is
carried out via the simplified formula
<\sigma^2> = 1/N trace(\Omega_U)
derived in our paper, where \Omega_U is the covariance matrix between the nodes
in the space orthogonal to the fully synchronized state vector \psi_0.
This is also equal to U^T \Omega U, the projection of the covariance matrix
\Omega via the unaveraging operator U, where U = I - G, G_ij = 1/N.

As such, there are two steps involved here:
1. Calculating the projected covariance matrix, and
1. Computing <\sigma^2> from that.

```matlab
% First compute the projected covariance matrix. discreteTime is a boolean
%  for whether we consider discrete-time (true) or continuous-time process,
%  MaxK is the number of iterations of the power series to compute \Omega_U
[sortedLambdasCU, Omega_U, err] = covarianceUGaussianNet(C, discreteTime, MaxK);
% Then compute sigma^2:
exp_sigma_sqr = synchronizability(Omega_U);
```

## Batch experiments involving parameter sweeps and repeat runs

In order to experiment with sweeping over many different parameters,
with many different sample networks, you can use (or build on) the
`computeSyncResults.m` script, which incorporates the above use cases
into an experimental framework.

First, set up a parameters file, using `parametersTemplate.m` as a template.
This will tell the compute script, for example, what size of network,
what type of network, the weights to use on coupled nodes, how many
repeat runs or samples for each parameter set to use, etc.
The sample parameters template is currently set up to recreate the results
for Fig 4b of the paper (i.e. N=100, d=4, c=0.5, p=[0.001, 0.002, 0.005, 0.01,
0.02, 0.05, 0.1, 0.2, 0.5, 1], directed network, continuous time,
no empirical runs (S=0), motif approximations up to 50.), except for only 2 repeat
runs per set (to make your first run go fast!).

Then call `computeSyncResults.m`, passing in the parameters file you have set up
(this can be supplied either as a filename, or running the script first
and then passing in `parameters` object it creates):

```matlab
computeSyncResults('parametersTemplate.m'); % using our template file
```

This call will run experiments across all of the specified parameters, 
with multiple network samples for each parameter set, and then
save the experimental results into a `.mat` file, in the 
folder specified by `parameters.folder`, formatted with a filename composed
from the parameters themselves. Running the above from our parametersTemplate.m
file will create a results file named `N100-randRing-d4-b1.00-c0.50-dir-k50-cont-S0-repeats2.mat`.

We show how to generate such batch results recreating the plots from our paper
in the next section.

## Running on a cluster

For small number of repeats or no empirical simulations, the scripts run 
fairly fast. To investigate large numbers or repeats or empirical simulations, particularly for large
numbers of samples, you'll likely be better served running the simulations
on a compute cluster, splitting out repeat runs over multiple compute processes.

See folder [cluster](tree/master/cluster) for suggestions on how to do this.

## Making plots from results files

After having generated results files, you can post-process these results to make plots.

To plot <\sigma^2> versus the swept parameter (p or c), you can call
`plotSyncResults.m`, passing in the same parameters object / filename:

```matlab
plotSyncResults('parametersTemplate.m'); % using our template file
```

Alternatively, plot the error between emprical and analytic calculations
for <\sigma^2> by calling `plotErrorInEmpiricalSyncResults.m`:

```matlab
`plotErrorInEmpiricalSyncResults('parametersTemplate.m'); % using our template file
```

This specifically requires having set the parameter S (number of empirical
samples) to be > 0 in order to have run empirical simulations, and having
run `computeSyncResults` one time for each of the values of `S` that you
wish to plot for (specified by `parameters.SRangeToPlot`).

We show how to generate such batch results recreating the Figure 1 and 4 plots
from the paper below.

# Recreating the results from our paper

The folder [2023-AnalyticRelationshipPaper](tree/master/2023-AnalyticRelationshipPaper)
contains scripts / parameters files to recreate these results.

## Figure 1 - convergence of numerical results from simulations

The parameter settings for recreating figure 1 are contained in `figure1Parameters.m`.

The script to run the experiments is `figure1ConvergenceOfSimulations.m`.
To run this, start in Matlab from the main toolkit folder and run:
```matlab
cd 2023-AnalyticRelationshipPaper
figure1ConvergenceOfSimulations
```

An `.eps` and `.mat` file is saved for the figure.

*Crucially* - note that the script adjusts the number of repeat runs 
downwards from 2000 in the parameters file (matching the paper) to 10 -- this
is done to ensure that the results run quickly for you. 
Also, the script removes the largest number of time series samples (L=1,000,000)
for the same reason.
These shorter runs take about 8 minutes 
on my machine -- running for 2000 repeats and with more time series samples
will take about 100 hours (the longest number of time series samples really 
blows this out)
so you would be advised to run for more repeats in parallel fashion as per the cluster
instructions below.

In any case, the shorter runs give a good idea of the main trends already,
though there are some fluctuations / std error differences compared to the longer runs.

## Figure 4 - small world sweeps

The parameter settings for recreating figure 4a are contained in `figure4aParameters.m`,
and these parameters are then adjusted by the experimental script for the other
sub-figures.

The script to run the experiments is `figure4SyncThroughSmallWorld.m`.
To run this, start in Matlab from the main toolkit folder and run:
```matlab
cd 2023-AnalyticRelationshipPaper
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
on my machine -- running for 1000 repeats will take about 20 hours,
so you would be advised to run for more repeats in parallel fashion as per the cluster
instructions below.

In any case, running 10 repeats only gives a good idea of the main trends already,
though there are some fluctuations / std error differences compared to the longer runs.

# Brief descriptions of files in this folder:

Primary experimental scripts to run:
* `computeSyncResults.m` - main script to sample <\sigma^2> for many networks
   with given parameters; saves results to .mat files for later processing.

Plotting scripts once results are ready:
* `plotSyncResults.m` - plot <\sigma^2> versus p or c (like Figure 4).
* `plotErrorInEmpiricalSyncResults.m` - plot difference between analytic
   and empirical results for <\sigma^2> (like Figure 1).

User-level scripts for analytical computation of projected covariance matrices and <\sigma^2>: 
* `covarianceUGaussianNet.m` - computes the projected covariance matrix (U^T \Omega U) and eigenvalues of a given connectivity matrix C.
* `synchronizability.m` - computes <\sigma^2> from the output (U^T \Omega U)  of `covarianceUGaussianNet.m`.

Underlying scripts involved in _analytical_ computation of projected covariance matrices and <\sigma^2>:
* `contCon2CovProjected.m` * - computes the projected covariance matrix for the continuous-time case
* `discreteCon2CovProjected.m` * - computes the projected covariance matrix for the discrete-time case

Underlying scripts to run _numerical_ simulations of dynamics for empirical calculation of <\sigma^2>:
* `empiricalCovariancesProjected.m` - runs numerical simulations to return an empirical measurement of the projected covariance matrix, for a given number of time samples S.
* `contCon2LaggedCovProjected.m` * - computes covariance between nodes, starting from zero covariance at time 0; used in the simulations of continuous time dynamics.

Scripts to generate new adjacency matrices:
* `generateNewRandomMatrix.m` - generate new Erdos-Renyi graph
* `generateNewRandomMatrixFixedD.m` - generate network with a fixed in-degree for each node, sources selected at random
* `generateNewRandomRingMatrix.m` - generate new Watts-Strogatz ring network, starting with fixed in-degree and then randomising edges (for directed graphs the in-degree is maintained)

Sample parameters file:
* `parametersTemplate.m` - set parameters for batch experimental runs.

Utility files:
* `adjMatrixToList.m` - converts an adjacency matrix to a list, required for some of the network generators.
* `parseParameters.m` - used by experimental scripts to convert varargin to pull out all required parameters.
* `maxrelerr.m` *  - used in convergence checks for the power series calculation of the projected covariance matrix
* `negligible.m` * - used in convergence checks for the power series calculation of the projected covariance matrix
* `runDfs.m` - runs a depth first search, used to ensure generated networks are (undirectionally) connected.

Files marked * are adapted from the ncomp toolkit of Barnett et al., 
originally distributed at http://www.secse.net/ncomp/ncomp_tools.zip
(though no longer available at that address), with the paper
L. Barnett, C. L. Buckley and S. Bullock (2009),
"On Neural Complexity and Structural Connectivity",
Physical Review E, 79, 051914,
under a license to "use as you wish".

# News

_3/3/2023_ - Initial code version uploaded (repo currently private)

# Acknowledgements

JL was supported through the Australian Research Council DECRA grant DE160100630, and The University of Sydney SOAR award.
