function runComputeSyncResults(parameters)
%
% Compute the sync results for randomly generated networks as specified in the 
% parameters object (or file).
% This script is mainly intended for running on a cluster, where
% we need to set the Matlab path also and then run over multiple values of
% S (number of empirical sample steps) set from parameters.SRangeToPlot
%
% Inputs:
% - parameters - an object containing the expected properties, or a string
%    describing the filename to run load this object in
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

if ischar(parameters)
    % Assume that this string contains a filename which when run will load
    % a properties object for this run
    eval(['run ', parameters]);
end
% Postcondition: parameters are in the parameters object

% Add the toolkit to the path
addpath(genpath(parameters.syncToolkitPath));

fprintf('Beginning calculation of sync results for folder %s\n', parameters.folder);

for S = parameters.SRangeToPlot
    tic
    parameters.S = S;
    computeSyncResults(parameters);
    toc
end


end
