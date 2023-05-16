# Linear Sync Toolkit (linsync)

Copyright (C) 2012- [Joseph T. Lizier](http://lizier.me/joseph/)

The `linsync` toolkit provides Matlab tools for analysing **synchronization** in networks of linearly coupled nodes.
Specifically, it provides implementations of the maths to measure and explore
the expected mean square deviation from synchronization $\left\langle \sigma^2 \right\rangle$ as a function of network coupling matrix $C$, with
the mathematical details provided in the following paper:

J.T. Lizier, F. Bauer, F.M. Atay, and J. Jost,
_"Analytic relationship of relative synchronizability to network structure and motifs"_,
Under submission,
2023

The tools compute $\left\langle \sigma^2 \right\rangle$ in both:
* Continuous-time Ornstein-Uhlenbeck dynamics: $$dX(t) = -X(t)(I-C)\theta dt + \zeta dw(t),$$ where
  $C$ is the update matrix, $w(t)$ is a multivariate Wiener process (uncorrelated) with covariance matrix $I$,
  $\theta$ is the reversion rate and $\zeta^2$ is noise strength, and
* Discrete-time auto-regressive processes: $$X(t+1) = X(t)C + \zeta R(t),$$ where
  where $C$ is the update matrix, $R(t)$ is uncorrelated mean-zero unit-variance Gaussian noise
  and $\zeta^2$ is noise strength.

In both cases the $N \times N$ weighted connectivity matrix $C$ above is specified in row-vector form,
the mathematics handles the case where the fully synchronized state vector
$\psi_0=[1,1,...,1]$ is an eigenvector of $C$ with eigenvalue $\lambda_0 = 1$,
and $$\left\langle \sigma^2 \right\rangle = \lim_{t \rightarrow \infty} \left\langle 1/N \sum_i (x_i(t) - \bar{x}(t) )^2 \right\rangle $$.
The above paper gives the solutions in continuous time:
$$ \left\langle \sigma^2 \right\rangle = \frac{\zeta^2}{2\theta} \sum_{m=0}^{\infty} & \frac{2^{-m}}{N} \sum_{u=0}^{m} \binom{m}{u} \mathrm{trace}\left( U (C^u)^T C^{m-u} U \right) $$,
and in discrete time:
$$ \left\langle \sigma^2 \right\rangle & = \frac{\zeta^2}{N} \sum_{u=0}^{\infty}{ \mathrm{trace}\left( U (C^u)^T C^{u} U \right)} $$.

Note that the toolkit fixes $\theta=\zeta=1$ since they act as constant multipliers in the above solutions;
if you want an answer for arbitray $\theta,\zeta$ then you should multiply the solution the software
returns for $\left\langle \sigma^2 \right\rangle$ as above.

Please **cite** your use of the toolkit via the above paper.

In this readme, we describe:
1. the [use cases](#1-use-cases) of the code and how to thread them together
1. how to [recreate](#2-recreating-the-results-from-our-papers) results from our papers
1. the [contents](#3-brief-descriptions-of-files-in-this-folder) of this top-level folder (the primary analysis scripts).

# 1. Use cases

In this section we briefly outline the primary use cases here, being:

1. [Generating sample network structures](#11-generating-network-structure) $C$ to investigate, and

2. [Computing the deviation from sync](#12-analysing-the-deviation-from-sync-for-a-given-network-c) $\left\langle \sigma^2 \right\rangle$ for a given network structure $C$.

Then building on those we discuss more involved use cases for:

3. [Batch experiments](#13-batch-experiments-involving-parameter-sweeps-and-repeat-runs) involving parameter sweeps and repeat runs, or

4. [Running experiments on a cluster](#14-running-on-a-cluster).

Finally, with the results generated you can:

5. [plot the results](#15-making-plots-from-results-files) from these runs.

## 1.1 Generating network structure

In order to compute $\left\langle \sigma^2 \right\rangle$ for a network, we need a directed weighted connectivity matrix $C$ for it.
Several scripts are distributed for the user to generate $C$ matrices for standard structures, including:
* Erdos-Renyi random graphs -- `generateNewRandomMatrix.m`
* Fixed in-degree random graphs -- `generateNewRandomMatrixFixedD.m`
* Watts-Strogatz ring networks with fixed in-degree -- `generateNewRandomRingMatrix.m`

These scripts generate unweighted adjacency matrices, then the user should weight the edges.
A sample generating a weighted, directed Watts-Strogatz ring network of $N$ nodes, with in-degree $d$,
and rewiring probability $p$, without self-connections allowed from re-wiring, and ensuring 
the network remains at least weakly-connected:

```matlab
% First generate the unweighted adjacency matrix without self-edges:
A = generateNewRandomRingMatrix(N, d, p, false, true, false, true);
% Then compute the in-degrees (should be all d in this example):
D = diag(sum(A));
% and weight the edges as c/d each and insert self-edges with weight b - c:
% (note that using b=1 ensures \psi_0 is an eigenvector of C with \lambda_0 = 1)
C = (b - c) .* I + c .* A * inv(D);
```

## 1.2 Analysing the deviation from sync for a given network _C_

Computing the expected deviation $\left\langle \sigma^2 \right\rangle$ from the synchronized state for
a given network connectivty matrix $C$ is
carried out via the simplified formula
$$\left\langle \sigma^2 \right\rangle = \frac{1}{N} \mathrm{trace}(\Omega_U)$$
derived in our paper, where $\Omega_U$ is the covariance matrix between the nodes
in the space orthogonal to the fully synchronized state vector $\psi_0$.
This is also equal to $U^T \Omega U$, the projection of the covariance matrix
$\Omega$ via the unaveraging operator $U$, where $U = I - G$, $G_{ij} = 1/N$.

As such, there are two steps involved here:
1. Calculating the projected covariance matrix $\Omega_U$, and
1. Computing $\left\langle \sigma^2 \right\rangle$ from that.

```matlab
% First compute the projected covariance matrix. The two latter arguments are:
%  discreteTime is a boolean for whether we consider discrete-time (true) or continuous-time process (false).
%  MaxK is the number of iterations of the power series to compute \Omega_U.
[sortedLambdasCU, Omega_U, err] = covarianceUGaussianNet(C, false, 100000);
% Then compute sigma^2:
exp_sigma_sqr = synchronizability(Omega_U);
```
The call to `covarianceUGaussianNet` will throw an error if the matrix does not meet synchronisation conditions,
or will not ensure guaranteed convergence of $\Omega_U$.

## 1.3 Batch experiments involving parameter sweeps and repeat runs

In order to experiment with sweeping over many different parameters,
with many different sample networks, you can use (or build on) the
`computeSyncResults.m` script, which incorporates the above use cases
into an experimental framework.

First, set up a parameters file, using `parametersTemplate.m` as a template.
This will tell the compute script, for example, what size of network,
what type of network, the weights to use on coupled nodes, how many
repeat runs or samples for each parameter set to use, etc.
The sample parameters template is currently set up to recreate the results
for Fig 4b of the paper (i.e. `N=100`, `d=4`, `c=0.5`, `p=[0.001, 0.002, 0.005, 0.01,
0.02, 0.05, 0.1, 0.2, 0.5, 1]`, directed network, continuous time,
no empirical runs (`S=0`), motif approximations up to 50.), except for only 2 repeat
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

We show how to generate such batch results recreating the plots from our 2023 paper
in [2023-AnalyticRelationshipPaper](/2023-AnalyticRelationshipPaper).

## 1.4 Running on a cluster

For small number of repeats or no empirical simulations, the scripts run 
fairly fast. To investigate large numbers or repeats or empirical simulations, particularly for large
numbers of samples, you'll likely be better served running the simulations
on a compute cluster, splitting out repeat runs over multiple compute processes.

See folder [cluster](/cluster) for suggestions on how to do this.

## 1.5 Making plots from results files

After having generated results files, you can post-process these results to make plots.

To plot $\left\langle \sigma^2 \right\rangle$ versus the swept parameter ($p$ or $c$), you can call
`plotSyncResults.m`, passing in the same parameters object / filename:

```matlab
plotSyncResults('parametersTemplate.m'); % using our template file
```

Alternatively, plot the error between emprical and analytic calculations
for $\left\langle \sigma^2 \right\rangle$ by calling `plotErrorInEmpiricalSyncResults.m`:

```matlab
plotErrorInEmpiricalSyncResults('parametersTemplate.m'); % using our template file
```

This specifically requires having set the parameter `S` (number of empirical
samples) to be $> 0$ in order to have run empirical simulations, and having
run `computeSyncResults` one time for each of the values of `S` that you
wish to plot for (specified by `parameters.SRangeToPlot`).

We show how to generate such batch results on a cluster recreating the Figure 1 and 4 plots
from the 2023 paper in [cluster](/cluster).

# 2. Recreating the results from our papers

* For Lizier et al., "Analytic relationship of relative synchronizability to
network structure and motifs", 2023, the folder [2023-AnalyticRelationshipPaper](/2023-AnalyticRelationshipPaper)
contains scripts / parameters files to recreate these results, as well as a README documenting how to run them.

# 3. Brief descriptions of files in this folder:

Primary experimental scripts to run:
* `computeSyncResults.m` - main script to sample $\left\langle \sigma^2 \right\rangle$ for many networks
   with given parameters; saves results to .mat files for later processing.

Plotting scripts once results are ready:
* `plotSyncResults.m` - plot $\left\langle \sigma^2 \right\rangle$ versus $p$ or $c$ (like Figure 4).
* `plotErrorInEmpiricalSyncResults.m` - plot difference between analytic
   and empirical results for $\left\langle \sigma^2 \right\rangle$ (like Figure 1).

User-level scripts for analytical computation of projected covariance matrices and $\left\langle \sigma^2 \right\rangle$: 
* `covarianceUGaussianNet.m` - computes the projected covariance matrix ($\Omega_U$) and eigenvalues of a given connectivity matrix $C$.
* `synchronizability.m` - computes $\left\langle \sigma^2 \right\rangle$ from the output ($\Omega_U$)  of `covarianceUGaussianNet.m`.

Underlying scripts involved in _analytical_ computation of projected covariance matrices and $\left\langle \sigma^2 \right\rangle$:
* `contCon2CovProjected.m` * - computes the projected covariance matrix for the continuous-time case
* `discreteCon2CovProjected.m` * - computes the projected covariance matrix for the discrete-time case

Underlying scripts to run _numerical_ simulations of dynamics for empirical calculation of $\left\langle \sigma^2 \right\rangle$:
* `empiricalCovariancesProjected.m` - runs numerical simulations to return an empirical measurement of the projected covariance matrix, for a given number of time samples `S`.
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
_"On Neural Complexity and Structural Connectivity"_,
Physical Review E, 79, 051914,
under a license to "use as you wish".

# 4. News

_3/3/2023_ - Initial code version uploaded (repo currently private)

# 5. Acknowledgements

JL was supported through the Australian Research Council DECRA grant DE160100630, and The University of Sydney SOAR award.

As above, several files were adapted from the ncomp toolkit of Barnett et al.
