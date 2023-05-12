function result = isdiagonalizable(C)
%% function result = isdiagonalizable(C)
%
% Quick check of whether the connectivity matrix is diagonlizable
%
% Inputs:
% - C - connectivity matrix to check
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

    tolerance = 1e-5;

    % Check that we can recompose directly from the left eigenvectors and eigenvalues (slower):
    [~, D, leftEigsTransposed] = eig(C);
    recomposedFromLeftEigs = inv(leftEigsTransposed') * D * leftEigsTransposed';
    if (max(max(abs(C - recomposedFromLeftEigs))) > tolerance)
        % The diagonalization didn't work well
        result = false;
    else
        result = true;
    end

end