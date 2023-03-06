function [empiricalCovarianceProjected] = empiricalCovariancesProjected(C, discreteTime, S, dt)
%
% Returns an empirical measurement of the covariance in the space orthogonal to the zero mode psi_0=[1 1 1 ... 1]
%  for either continuous time Ornstein-Uhelnbeck or discrete time AR
%  dynamics as required.
% 
% There was no original implementation of this in the Barnett ncomp
% toolkit, so we built the Cholesky decomposition here following their
% guidelines. Differently to that, we are projecting orthogonal to the
% fully synchronized mode at every time step, since we want the projected
% covariance matrix in the end anyway and this will give a more stable
% result.
%
% Inputs:
% - C - network coupling matrix (C(i,j) = i->j weight, in row vector form)
% - discreteTime - whether the time-series is estimated from a discretized process (true) or use of exact method assuming continuous process (false).
% - S - number of time steps (sample size) of the time series to compute the covariance over.
% - dt - time interval for iterating the exact method (continuous time). Ignored for discrete time.
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

N = size(C,1);
N2 = size(C,2);
if (N2 ~= N) 
    error('C is not square');
end
I = eye(N);
% Averaging projection operator:
G = 1 ./ N * ones(N);
% Unaveraging operator:
U = I - G;

% Burn-in time before keeping samples for analysis:
initialTimeSteps = S; % For no lambda_0=1, this parameter doesn't seem to matter much, it settles quickly. No difference using 1000 or 100000 steps. For having \lambda_0=1, having at least 10000 steps would be advisable, though the accuracy seems to depend more on S than this.

tol = 10; % Multiples of machine epsilon within which we want to consider a value to be zero.

if (discreteTime)
    % Run the discrete time multivariate AR (autoregressive) process

    X = normrnd(0, 1, 1, N); % initialise row vector, with std deviation of the noise
    Xproj = X * U; % Projection away from the synchronized state
    allXproj = zeros(S,N); % Samples of X * U
    sampleNum = 1;
    for t = 1 : initialTimeSteps + S
        % Make the next time step:
        randNoise = normrnd(0, 1, 1, N);
        % Since we're going to project by U in the end anyway to compute sync, we're going
        % to make the projection at every step of the dynamics - this
        % avoids the values of X alone becoming unstable and compromising
        % the final empirical calculation.
        newXproj = Xproj * C * U + randNoise * U;
        if (t > initialTimeSteps)
            % Keep the new values
            allXproj(sampleNum, :) = newXproj;
            sampleNum = sampleNum + 1;
        end
        Xproj = newXproj;
    end
else
    % Run the Ornstein-Uhlenbeck process using exact method
    
    % How do the expected values of X adjust over dt: (appendix A of Barnett et al)
    %  <X(t+s)> = X(t) * meanMultiplier
    meanMultiplier = expm(-(I-C)*dt);

    % Generate the covariance matrix M(dt) for how the variables X covary
    %   after we make the dt time step if they started with no covariance:
    symmetricC = isempty(find(C - C' > tol*eps));
    [MsProjected,err] = contCon2LaggedCovProjected(C,dt,symmetricC,1000,tol,true);
    if (err ~= 0)
        return;
    end
    % Make a Cholesky decomposition - returns the upper triangle, or if MsProjected is rank deficient,
    %  which it is because we project orthogonal to the zero mode psi_0, then (from matlab help):
    % "R is an upper triangular matrix of size q-by-n so that the L-shaped region of the first q rows and first q
    %  columns of R'*R agree with those of A."
    [cholR,p] = chol(MsProjected);
    fprintf('chol returned p=%d, giving q=%d rows and cols in cholR\n', p, p - 1);
    if (p > 0)
        % Pad the cholesky decomposition with zeros to make an NxN matrix to multiply the noise by.
        % We assume that the zero multiplications here still have something come through for those columns
        %  because of the unaveraging procedure.
        noiseMultiplier = [cholR, zeros(p - 1, N - (p - 1)); zeros(N - (p - 1), N)];
    else
        noiseMultiplier = cholR;
    end

    X = normrnd(0, 1, 1, N) * sqrt(dt); % row vector like Barnett and Bullock, with std dev of the noise to start with
    Xproj = X * U;
    allXproj = zeros(S,N);
    sampleNum = 1;
    for t = 1 : initialTimeSteps + S
        % Make the next time step:
        % Since we're going to project by U in the end anyway to compute sync, we're going
        % to make the projection at every step of the dynamics - this
        % avoids the values of X alone becoming unstable and compromising
        % the final empirical calculation.

        % The mean of each vector after projection will be:
        newXProjMean = Xproj * meanMultiplier * U;
        % and the updates to the process over s will have covariance matrix MsProjected,
        %  so we work out the variances via a Cholesky decomposition (App A of Barnett):
        newXProj = newXProjMean + normrnd(0, 1, 1, N) * noiseMultiplier * U;
        if (t > initialTimeSteps)
            % Keep the new values
            % allX(:, t - multiplesOfSForSimulating*S) = newX; % If we're doing the S after 10*S
            allXproj(sampleNum, :) = newXProj;
            sampleNum = sampleNum + 1;
        end
        Xproj = newXProj;
    end
end

% If we had not made the projection at every step of the dynamics above, we would now as follows:
% empiricalCovarianceProjectedFromCov = cov(allX * U);
% But in theory the above may diverge (though in practise we won't be using enough samples for that to happen)
%  because the time series may have diverged before we can take the
%  covariance. So instead we project at every time step above, achieves
%  same as projecting now but is more stable.

empiricalCovarianceProjected = cov(allXproj);

end

