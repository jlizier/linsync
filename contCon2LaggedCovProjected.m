function [M_U,err] = contCon2LaggedCovProjected(C,s,symmetric,maxi,tol,verbose)
%
% This code computes the projected lagged covariance matrix for the contimous-time
% Ornstein-Uhlenbeck case after time s, assuming starting from no
% correlation at time 0.
% The code converts con2lcov from ncomp from Barnett et al to handle
%  the projection of the covariance matrix orthogonal to the synchronized
%  vector.
% See reference for the Barnett code:
%   - ncomp_tools: L. Barnett, C. L. Buckley and S. Bullock (2009)
%      (licensed to "use as you wish")
%   - L. Barnett, C. L. Buckley and S. Bullock (2009)
%       On Neural Complexity and Structural Connectivity, Physical Review E, 79, 051914
%
% Calculates the lagged covariance matrix (after no covariance), projected into the space orthogonal to the
%  zero mode eigenvector \psi_0 = [1,1,1,1,...,1] for the stationary multivariate
%  Ornstein-Uhlenbeck process:
%
%     dX(t) = -X(t)*(I-C)s + dW_t
%
% where W_t is a multivariate Wiener process with identity covariance
% matrix and C is the connection matrix (1st index efferent)
%
% The notes here call this the lagged covariance, but this doesn't mean cov(x(t),x(t+s)) across a lag,
% it means M(s) = cov(x(t+s),x(t+s)) with s the delay since there was zero covariance (M(0) = 0).
% M_U is the projected version of M(s): M_U = M(s) * U
%
% The code solves the matrix ODE: dM_U/ds = U + U*C'*M_U + M_U*C*U - 2*M_U
% with inital condition M(0) = 0 and where M_U = U * M * U
% This adapts the ncomp code to project by U orthogonal to the fully
% synchronized mode.
% See Appendix A of the Barnett paper above,
%
% Inputs
% - C - connectivity matrix, including self-connection weights (in the Cij = i -> j format, i.e. for row-vectors)
% - s - is the lag (needn't be small)
% - symmetric - if C is symmetric to use (faster) non-iterative formula.
%      Currently this is not implemented so the argument is ignored.
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
% - verbose - if true, the code reports the number of iterations and convergence results
%
% Outputs
% - M_U - covariance matrix of the process after s steps projected into the orthonormal
%    space to the fully synchronized state vector psi_0
% - err - 0: no error, 1: warning that maximum iterations were reached, 2: covergence failed
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

err = 0;
n = size(C,1);
I = eye(n);
% Averaging projection operator:
G = 1 ./ n * ones(n);
% Unaveraging operator:
U = I - G;

% TODO Could return a proper symmetric form -- Barnett code was incomplete,
%  some work towards it below, but we've never really needed it:
%if symmetric % explicit solution available
%    % JL the original line here was:
%    % M  = (I-expm(-s*H))/(2*(I-C));
%    % BUT THIS DIDN'T DEFINE H!!
%    % So we correct it to match eqn A2 in their paper:
%    M_U  = (I-expm(-2*s*(I-C)))/(2*(I-C));
%    return
%end

M_U = U*s;
dM_U = M_U;
UCprime = U*C';
CU = C*U;
for i = 2:maxi
    dM_U = (UCprime*dM_U+dM_U*CU-2*dM_U)*(s/i);
    if negligible(dM_U,M_U,tol)
       break
    end
    M_U = M_U + dM_U;
end

mre = maxrelerr(dM_U,M_U);

if isnan(mre) || isinf(mre)
    fprintf(2,'ERROR (contCon2LaggedCovProjected): failed to converge\n');
    err = 2;
    return;
end

if i == maxi
    mac = mean(mean(abs(M_U)));
    fprintf(1,'WARNING (contCon2LaggedCovProjected): exceeded maximum iterations: mean abs. covariance = %d, max. relative error = %d\n',mac,mre);
    err = 1;
    return;
end

if verbose
    fprintf(1,'contCon2LaggedCovProjected: iterations = %i, max. relative error = %d\n',i,mre);
end

