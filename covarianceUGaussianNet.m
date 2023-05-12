function [sortedLambdasCU, UcovarianceU, err] = covarianceUGaussianNet(C, discreteTime, maximumIterations, forcePowerSeriesForSymmetric, verbose, skipPowerSeriesConvergenceCheck)
%
% Computes the eigenvalues of C, and the projected covariance matrix (U^T \Omega U) for the given network.
%
% Inputs
% - C - connectivity matrix, including self-connection weights (in the Cij = i -> j format, i.e. for row-vectors)
% - discreteTime - whether the time-series is estimated from a discrete-time AR process (true)
%    or use of exact method assuming continuous Ornstein-Uhlenbeck process.
% - maximumIterations - the number of components to add into the power series for UcovarianceU. Default is 1000
% - forcePowerSeriesForSymmetric - force the use of power series for UcovarianceU even if the
%     connectivity matrix is symmetric
% - verbose - level of verbosity for discreteCon2CovProjected and contCon2CovProjected
%     (set to -1 if only using a small number of iterations, to avoid
%     warnings)
% - skipPowerSeriesConvergenceCheck - to skip checking that max eigenvalue
%     lies within unit circle (only has an impact for continuous time
%     dynamics). Default is false.
% 
% Outputs
% - sortedLambdasCU - sorted eigenvalues of CU (from smallest to largest magnitude if complex)
% - UcovarianceU - projected covariance matrix U^T \Omega U, for U = I - 1/N and
%     \Omega is the covariance matrix, assuming the Ornstein-Uhlenbeck or AR process.
% - err - error status from the call to work out the covariance matrix here
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3


n = size(C,1);
I = eye(n);
G = ones(n) / n;
U = I - G;

if (nargin < 3)
    maximumIterations = 1000;
end
tol = 100; % Multiples of machine epsilon within which we want to consider a value to be zero.

if (nargin < 4)
    forcePowerSeriesForSymmetric = false;
end

if (nargin < 5)
    verbose = 0;
end

if (nargin < 6)
    skipPowerSeriesConvergenceCheck = false;
end

% Compute eigenvalues of C * U:
%  (since C is square, the eigenvalues are the same as C^T - i.e. it doesn't matter that they correspond to row/column vectors)
lambdasCU = eig(C * U);
sortedLambdasCU = sort(lambdasCU); % Sorts the elements by magnitude
% Find the minimum and maximum eigenvalues
lambdaCUMax = sortedLambdasCU(size(sortedLambdasCU,1));
    
% Check for stationarity for C*U
if (discreteTime)
    if (abs(lambdaCUMax) >= 1)
        % Non-stationary
        save('nonconvergentNetwork.mat', '-mat', 'C'); % save for later investigation
        error('Discrete system with |\\lambda_CU_max| >= 1 (%.6f) i.e. non-stationary of C*U\n', abs(lambdaCUMax));
    end
else
    sortedRealPartsOfLambdasCU = sort(real(lambdasCU));
    realLambdaCUMax = sortedRealPartsOfLambdasCU(size(sortedRealPartsOfLambdasCU,1));
    if (realLambdaCUMax >= 1)
        % Non-stationary
        save('nonconvergentNetwork.mat', '-mat', 'C'); % save for later investigation
        error('Continuous system with Max(Re(\\lambda_CU)) >= 1 (%.6f), i.e. non-stationary of C*U\n', realLambdaCUMax);
    end
end

% Check for convergence of the projected covariance matrix: (same condition for both continuous and discrete):
if  (~skipPowerSeriesConvergenceCheck && (abs(lambdaCUMax) >= 1))
    save('nonconvergentNetwork.mat', 'C'); % save for later investigation
    error('|\\lambda_CU_max| >= 1 (%.2f) implies that the U^T \\Omega U matrix will not converge.\n', abs(lambdaCUMax));
end

% Check whether matrix C is symmetric or not, to help speed up the
% projected covariance calculation:
symmetric = isempty(find(C - C' > tol*eps));

% Now compute the UcovarianceU matrix
if (discreteTime)
    [UcovarianceU,err] = discreteCon2CovProjected(C,symmetric && ~forcePowerSeriesForSymmetric,maximumIterations,tol,verbose);
else
    [UcovarianceU,err] = contCon2CovProjected(C,symmetric && ~forcePowerSeriesForSymmetric,maximumIterations,tol,verbose);
end

% Can handle err == 2 here, or just let it go through to the caller

end

