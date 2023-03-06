function [UcovarianceU,err] = discreteCon2CovProjected(C,useSymmetricForm,maxi,tol,verbose)
%
% This code computes the projected covariance matrix for the discrete-time
% case.
% The code converts con2cov from ncomp from Barnett et al to the discrete case
%  and for the projection.
%  See reference for the Barnett code:
%   - ncomp_tools: L. Barnett, C. L. Buckley and S. Bullock (2009)
%      (licensed to "use as you wish")
%   - L. Barnett, C. L. Buckley and S. Bullock (2009)
%       On Neural Complexity and Structural Connectivity, Physical Review E, 79, 051914
%
% Calculates the projected covariance matrix U \Omega U for the multivariate
% discrete AR process (where X may have eigenvector \psi_0 but is otherwise stationary):
%
%     X(t+1) = X(t)*C + R(t)
%
% where R(t) is uncorrelated mean-zero unit-variance Gaussian noise
%  and C is the update matrix.
% The projection is orthogonal to the fully synchronized state vector psi_0.
%
% C(i,j) represents the update effect of i -> j ;
%  Note this follows the format of Barnett et al, where activity is in row vectors.
%
% Inputs
% - C - connectivity matrix, including self-connection weights (in the Cij = i -> j format, i.e. for row-vectors)
% - useSymmetricForm - whether C is symmetric or not AND whether we should use the
%       analytic shortcut for symmetric C's. Should set this to false if
%       you want to see apparoximations to a limited order.
% - maxi - maximum number of iterations.
% - tol - tolerance for concluding convergence has been reached.
%  These parameters are used the same way as it the ncomp toolbox of Barnett et al:
% "
%     'maxi' gives the maximum iterations for the non-symmetric version.
%     Convergence criterion is to stop when the ratios of the added terms
%     approach the machine epsilon with given tolerance 'tol'. If maximum
%     iterations are exceeded a warning is issued. Failure to converge
%     is treated as an error. You can play off accuracy against no. of iterations.
%     Try maxi = 1000, tol = 10 for starters. If it converges you should get a max.
%     relative error of order 1e-15 or better.
% "
% - verbose - whether to report warnings and the number of iterations and convergence results
%   - -1 - minimal, output errors only, not warnings or number of iterations and convergence results
%        (useful for batch runs with low order approximations)
%   - 0 (or false) - default, output errors and warnings but not number of iterations and convergence results
%   - 1 - verbose, output errors, warnings, and number of iterations and convergence results
%
% Outputs
% - UcovarianceU - covariance matrix of the process projected into the orthonormal
%    space to the fully synchronized state vector psi_0
% - err - 0: no error, 1: warning that maximum iterations were reached, 2: covergence failed
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

err = 0;
N = size(C,1);
I = eye(N);
% Averaging projection operator:
G = 1 ./ N * ones(N);
% Unaveraging operator:
U = I - G;

if (useSymmetricForm)
    % Simplest expression would be: UcovarianceU = U * inv(I - C*C); HOWEVER:
    %  this only holds where max(lambda) != 1.
    % A more general expression which does work for max(lambda) = 1 is:
    %  (which works because the denominator does not have the lambda=1 
    %  eigenvalue anymore) -- see paper appendix C.B
    UcovarianceU = inv(I - (C*U)^2)*U;
    return;
end

UcovarianceU = U; % Technically this is transpose(U) * U, but this is the same as U
leftComponent = U;
rightComponent = U;
leftMultiplier = U * C';
rightMultiplier = C * U;

for i = 1:maxi
    % dProjM holds the previous C'^i * C^i term
    % ProjM holds the sum of previous projected C'^i * C^i terms
    leftComponent = leftComponent * leftMultiplier;
    rightComponent = rightMultiplier * rightComponent;
    d_UcovarianceU = leftComponent * rightComponent;
    if negligible(d_UcovarianceU,UcovarianceU,tol)
       break
    end
    UcovarianceU = UcovarianceU + d_UcovarianceU;
end

mre = maxrelerr(d_UcovarianceU,UcovarianceU);

if isnan(mre) || isinf(mre)
    fprintf(2,'ERROR (discreteCon2CovProjected): failed to converge\n');
    err = 2;
    return;
end

if (i == maxi) && (verbose >= 0)
    mac = mean(mean(abs(UcovarianceU)));
    fprintf(1,'WARNING (discreteCon2CovProjected): exceeded maximum iterations (%d): mean abs. covariance = %.3f, max. relative error = %.3e (both in projected space)\n',maxi,mac,mre);
    err = 1;
    return;
end

if (verbose > 0)
    fprintf(1,'discreteCon2CovProjected: iterations = %i, max. relative error = %d\n',i,mre);
end

