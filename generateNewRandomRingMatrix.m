function A = generateNewRandomRingMatrix(N, d, gamma, includeSelf, undirected, allowSelf, ensureConnected)
%% Generate a new ring network, which is rewired subject to a given probability
%
% Inputs
% - N - network size
% - d - in-degree
% - gamma - probability of rewiring each link after the regular network is constructed
% - includeSelf - whether to include self-connections in the regular
%    network (included in degree d)
% - undirected - whether to make the matrix directed or not.
%   For directed graphs, we rewire the sources, to keep a fixed in-degree
%   For undirected graphs, we don't allow connection to an odd number of
%   other nodes (too unconstrained to work out who to connect to)
% - allowSelf - whether to allow rewired connected to be made to oneself.
% - ensureConnected - only return an at least weakly-connected matrix.
%     If this is set to true, then p must be >= 2/N(N-1) for undirected, 1/N(N-1) for directed, for the matrix to be connected.
%     Note: if allowSelf==true, then (N-1) -> N here.
%
% Outputs
% - A - connectivity matrix (can be directed; A(i,j) means a link exists from i->j)
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

fprintf('Running\n');
if (includeSelf)
    d = d - 1;
end

while (true)
    % Generate a new regular network structure with degree d.
    % Start by assigning self-connections if required.
    if (includeSelf)
        A = eye(N);
    else
        A = zeros(N);
    end        

    % Now assign the d links for each node
    if (mod(d, 2) == 0)
        % we need to assign an even number of links. Just assign one
        % direction for the moment, and then we'll mirror them
        for l = 1 : d/2
            A = A + diag(ones(N-l,1), l) + diag(ones(l,1), -(N-l));
        end
        % Now mirror the undirected edges so that we have them on both
        % sides of each node
        A = A | A';
    else
        % We need to assign an odd number of links.
        if undirected
            error('We do not handle an odd number of edges for undirected ring networks\n');
        end
        % Postcondition - we've got a directed network
        % First assign one direction, and then we'll mirror them.
        for l = 1 : floor(d/2)
            A = A + diag(ones(N-l,1), l) + diag(ones(l,1), -(N-l));
        end
        % Now mirror the undirected edges for these inner links
        A = A | A';
        % And now add one extra incoming edge, which will come from the
        % lower index node (without loss of generality).
        % We only add an extra link in
        l = ceil(d/2);
        A = A + diag(ones(N-l,1), l) + diag(ones(l,1), -(N-l));        
    end
    % Error check: confirm that every node has the correct degree
    if  (includeSelf)
        if (find(sum(A) - ones(1,N)*(d+1)))
            error('Some nodes do not have degree %d\n', d+1);
        end
    else
        if (find(sum(A) - ones(1,N)*d))
            error('Some nodes do not have degree %d\n', d);
        end
    end

    % Now rewire the edges
    % First take a snapshot of which the original edges were:
    L = adjMatrixToList(A, true);
    edgesConsidered = 0;
    rewirings = 0;
    for dest = 1 : N
        for source = L{dest}
            if (undirected && (source < dest))
                % For undirected edges, try to only attempt to rewire them once
                % (it's possible an attempt could be made on a rewired link
                % if it gets moved to a node with higher index, so this isn't perfect)
                continue;
            end
            edgesConsidered = edgesConsidered + 1;
            if (rand < gamma)
                % Rewire this edge:
                rewirings = rewirings + 1;
                if (not(undirected))
                    % For directed graphs, rewire the sources (to keep fixed indegrees)
                    remainingNode = dest;
                    departingNode = source;
                else
                    % For undirected graphs, select the remaining node at random
                    if (rand < 0.5)
                        remainingNode = dest;
                        departingNode = source;
                    else
                        remainingNode = source;
                        departingNode = dest;
                    end
                end
                % Select the new candidate for connection into remainingNode
                candidates = find(A(:,remainingNode) == 0);
                if (not(allowSelf))
                    % Then make sure we can't become connected to ourself
                    selfCandidateIndex = find(candidates == remainingNode);
                    if (not(isempty(selfCandidateIndex)))
                        candidates(selfCandidateIndex) = [];
                    end
                end
                candidateIndex = ceil(rand * length(candidates));
                newNode = candidates(candidateIndex);
                % Now rewire the connection
                A(departingNode, remainingNode) = 0;
                A(newNode, remainingNode) = 1;
                if (undirected)
                    % Rewire the complementary directed edge also
                    A(remainingNode, departingNode) = 0;
                    A(remainingNode, newNode) = 1;
                end
                %printf('Rewired %d -> %d to %d -> %d\n', departingNode, remainingNode, newNode, remainingNode);
            end
        end
    end
    %printf('Rewired %d out of %d edges considered\n', rewirings, edgesConsidered);

    if (ensureConnected)
        % Now check that the matrix is at least weakly connected
        symmA = A | A';
        visited = zeros(1,N);
        visited = runDfs(symmA,N,1,visited);
        if (sum(visited) ~= N)
            % The matrix is not undirectionally connected
            printf('Matrix is not unidirectionally connected, trying again ...\n');
            continue;
        end
    end
    % The matrix is fine to return:
    % Debug print:
    actualIndegree = sum(sum(A)) ./ N;
    fprintf('Generated a ring graph with in-degree %.4f for selected d=%d, with %d rewirings out of %d edges considered\n', ...
        actualIndegree, d, rewirings, edgesConsidered);

    return;
end
end

