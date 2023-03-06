% Compute <sigma^2>, it's low order approximations, and leading eigenvalues
% for continuous-time Ornstein-Uhelnback process through a small-world
% transition.
%
% This recreates the plots in Figure 4 of our paper.
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

%% Preliminaries:

% Load in the parameters object for the first experiment, figure 4a:
figure4aParameters;

% Add the toolkit to the path:
addpath(genpath(parameters.syncToolkitPath));

% Make sure the folder for the results to be stored in exists:
if (exist(parameters.folder, 'dir') == 0)
    % Need to create the results folder:
    [success, message] = mkdir(parameters.folder);
    if (success == 0)
        error('Creating the results folder failed: %s\n', message);
    end
end

% To speed up the experiment when not on a cluster you can set:
parameters.repeats = 10; % instead of 1000
% You can get a decent idea of the main trends for this number already,
% though there are significant fluctuations / std error differences
% compared to the longer runs.

%% Figure 4a experiment and plots:
% Here, the very regular graphs (when there's at least one re-wiring to make
%  them asymmetric) take many iterations to converge, because of the large
%  connection strengths for only d=2 and large eigenvalues with the regular
%  graphs of these.
computeSyncResults(parameters);
% Takes about 2.5 mins on my machine for 10 repeats, 4 hours for 1000
% repeats. Better to do on cluster for that many.
% Now plot the results:
plotSyncResults(parameters);
print('-depsc', [parameters.folder, 'fig4a.eps'])
saveas(gca, [parameters.folder, 'fig4a.fig'], 'fig');

%% Figure 4b experiment and plots:
parameters.d = 4; % This is the only parameter that is different to 4a
% Converges faster than 4a because of smaller c/d weights.
computeSyncResults(parameters);
% Takes about 1.5 mins on my machine for 10 repeats, 2.5 hours for 1000
% repeats. Better to do on cluster for that many.
% Now plot the results:
plotSyncResults(parameters);
print('-depsc', [parameters.folder, 'fig4b.eps'])
saveas(gca, [parameters.folder, 'fig4b.fig'], 'fig');

%% Figure 4c experiment and plots:
parameters.d = 8; % This is the only parameter that is different to 4b
% Converges faster than 4a/b because of smaller c/d weights.
computeSyncResults(parameters);
% Takes about 0.75 mins on my machine for 10 repeats, 1.25 hours for 1000
% repeats. Better to do on cluster for that many.
% Now plot the results:
plotSyncResults(parameters);
print('-depsc', [parameters.folder, 'fig4c.eps'])
saveas(gca, [parameters.folder, 'fig4c.fig'], 'fig');

%% Figure 4d experiment and plots:
parameters.d = 4;
parameters.c = 0.1; % These are the only two parameters different to 4c
% Converges more slowly than 4a/b/c because of larger 1-c self-weights.
computeSyncResults(parameters);
% Takes about 5 mins on my machine for 10 repeats, 8.5 hours for 1000
% repeats. Better to do on cluster for that many.
% Now plot the results:
plotSyncResults(parameters);
print('-depsc', [parameters.folder, 'fig4d.eps'])
saveas(gca, [parameters.folder, 'fig4d.fig'], 'fig');

%% Figure 4e experiment and plots:
parameters.c = 1.0; % This is the only parameter different to 4d
% Converges faster than 4a-d because of no self-links and small c/d cross-weights.
computeSyncResults(parameters);
% Takes about 0.8 mins on my machine for 10 repeats, 1.3 hours for 1000
% repeats. Better to do on cluster for that many.
% Now plot the results:
plotSyncResults(parameters);
print('-depsc', [parameters.folder, 'fig4e.eps'])
saveas(gca, [parameters.folder, 'fig4e.fig'], 'fig');

%% Figure 4f experiment and plots:
% Switch to a sweep across c for fixed p:
parameters.p = 0.001;
parameters.c = 0.1:0.1:1.0;
computeSyncResults(parameters);
% Takes about 1.25 mins on my machine for 10 repeats, 2 hours for 1000
% repeats. Better to do on cluster for that many.
% Now plot the results:
plotSyncResults(parameters);
print('-depsc', [parameters.folder, 'fig4f.eps'])
saveas(gca, [parameters.folder, 'fig4f.fig'], 'fig');

fprintf('All experiments finished, using %d network samples (is this the same as the 1000 used for the paper?)\n', ...
    parameters.repeats);
