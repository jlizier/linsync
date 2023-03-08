%% Parses varargin into the relevant parameters
% either based on a structure, filename, or
% each parameter supplied individually.
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

if (length(varargin) == 1)
    % User has provided a parameters object directly or a string specifying
    % the filename to load a parameters object in.
    parameters = varargin{1};
    if ischar(parameters)
        % Assume that this string contains a filename which when run will load
        % a properties object for this run
        eval(['run ', parameters]);
    end
    % Postcondition: parameters are in the parameters object
    
    % Assign all of the relevant variables:
    N = parameters.N;
    b = parameters.b;
    c = parameters.c;
    undirected = parameters.undirected;
    discretized = parameters.discretized;
    repeats = parameters.repeats;
    networkType = parameters.networkType;
    p = parameters.p;
    d = parameters.d;
    S = parameters.S;
    SRangeToPlot = parameters.SRangeToPlot;
    maxMotifLength = parameters.maxMotifLength;
    folder = parameters.folder;
    MaxK = parameters.MaxK;
    dt = parameters.dt;
    randSeed = parameters.randSeed;
elseif (length(varargin) < 12)
    fprintf('Not enough arguments supplied, see code for details');
else
    % All parameters have been supplied individually
    N = varargin{1};
    b = varargin{2};
    c = varargin{3};
    undirected = varargin{4};
    discretized = varargin{5};
    repeats = varargin{6};
    networkType = varargin{7};
    p = varargin{8};
    d = varargin{9};
    S = varargin{10};
    SRangeToPlot = S; % only used by plot scripts in position of S
    maxMotifLength = varargin{11};
    folder = varargin{12};
    if (length(varargin) > 12)
        MaxK = varargin{13};
    else
        MaxK = 2000;
    end
    if (length(varargin) > 13)
        dt = varargin{14};
    else
        dt = 1;
    end
    if (length(varargin) > 14)
        randSeed = varargin{15};
    else
        randSeed = 'shuffle';
    end
end
