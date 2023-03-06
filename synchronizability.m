function exp_sigma_sqr = synchronizability(UcovarianceU)
% 
% Compute the expectation value of \sigma^2, meaning the expected average
% square deviation of nodes from the network mean.
%
% Inputs
% - UcovarianceU - covariance matrix, projected into space orthogonal to the
%    synchronized state psi_0, computed from the coupling matrix C
%
% Outputs
% - exp_sigma_sqr - expectation of the average square distance of nodes
%    from the network mean, a.k.a the width of the sync landscape
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

N = size(UcovarianceU,1);

exp_sigma_sqr = trace(UcovarianceU) ./ N;

end

