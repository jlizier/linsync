%% Script to set up the parameters object for the figure 4b experiment (d=4, c=0.5, various p).
% You can not only assign values here, but have differential processing
% (e.g. to do different things on your desktop or cluster).
% Required members are described as they appear below.
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

clear parameters; % In case parameters previously held the name of the parameters file

% You can use this boolean to flag different folders for local and cluster runs:
if exist('/home/joseph/') > 0
    isCluster = false;
else
    isCluster = true;
end

% Set the location of the code library - only required for
% runComputeSyncResults which is designed for cluster use
if isCluster
    parameters.syncToolkitPath = '..'; % Could be elsewhere on your cluster depending on how you've assigned folders
else
    parameters.syncToolkitPath = '..'; % Toolkit is in the parent folder
end
 
% - N - Network size
parameters.N = 100;

% - b - b and c are weighting factors for the variational equation from Atay et al, SIAM 2006.
%   b is the sum of all incoming weights for each node, including self.
%   So (b - c) is the strength of coupling to previous state of the destination.
%   For stability we require |b| < 1 when undirected, though b == 1 ensures
%   zero-mode as an eigenvector (for our first sync paper we only use b =
%   1).
parameters.b = 1.0;
% - c - as above; c is the total strength of non-self coupling of the sources to the destinations.
%   For stability we require |b - 2c| < 1 when undirected (though I don't
%   think this guarantees stability if undirected).
%   Can be an array to run experiments for each value (one of c or p can be
%   an array)
parameters.c = 0.5;

% - undirected - whether it is an undirected network or not.
parameters.undirected = false;

% - discretized - whether the time-series is estimated from a discrete-time
%    process (true) or use of exact method assuming continuous-time process (false)
parameters.discretized = false;

% - repeats - how many networks to sample with each parameter set.
%   Want to have this as a large number for good sampling.
%   Set it to a small number if you want test runs with this script to run
%   quickly, or on the cluster
parameters.repeats = 10;

% - networkType. Options are:
%   - 'rand' - random network
%   - 'randFixedD' - random network with fixed in-degree
%   - 'randRing' - a ring network with randomly rewired edges (Watts-Strogatz model).
parameters.networkType = 'randRing';

% Next 2 arguments that follow depend on which network type was requested:

% - p - is interpreted as either: connection probability for rand or randFixedD
%   (disconnected networks are rejected), or rewiring probability for randRing
% Can be an array (instead of c), such that we will run over all values of p.
parameters.p = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0];

% - d - degree for randRing network type. d/2 links on either side of the
%   node.
parameters.d = 4;

% - S - number of samples for the empirical calculation
%   (known as L in our first sync paper). Need very large
%   if you want high accuracy; see e.g. Fig 4 in Barnett (2009), though
%   with large \lambda's (e.g. for regular network) our errors will be
%   larger again. Try at least 1000.
%   If S==0, then we won't run any empirical calculations, just make the
%   analytic calculations.
%   Must be a scalar; put arrays for iterating over multiple values
%   into SRangeToPlot if using runComputeSyncResults etc
parameters.S = 0;
%   runComputeSyncResults (on cluster) or plotErrorInEmpiricalSyncResults
%   use the full range in the following:
parameters.SRangeToPlot = 0;


% - maxMotifLength - for continuous process, this is the max motif length (m)
%    to make the sync approximations up to; for discrete process, this is the maximum walk
%    length (u) to make the sync approximations up to (this is half the motif size,
%    since the relevant motif here is two converging walks of the same size u).
parameters.maxMotifLength = 50;

% - folder - directory where all of the files are to be stored.
if isCluster
    parameters.folder = './results/N100-randRing-d4-b1.00-c0.50-dir-k4-cont/[@P1]';
else
    parameters.folder = './results/';
end

% combineResultsFrom - in the case of a cluster run, we will want to combine
%  the results which are in various files in folders under this folder.
% They will be saved back to the parameters.folder above
parameters.combineResultsFrom = '~/temp/sync/empirical/';

% MaxK - maximum number of iterations for solving the power series for the
% covariance matrix.
% For the regular networks with small d (and therefore higher c/d), need up
% to tens of thousands of iterations; seems safe at 100000.
% Can leave it like this for all runs; they will terminate early when
% converging anyway.
parameters.MaxK = 100000;

% dt - time interval for iterating the exact method (continuous time). Ignored for discrete time.
parameters.dt = 1; 
