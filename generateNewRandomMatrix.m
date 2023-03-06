function A = generateNewRandomMatrix(N, p, allowSelf, undirected, ensureConnected, ensureALinkForEachNode)
%% 
% Generate a new Erdos-Renyi random matrix
%
% Inputs
% - N - network size
% - p - connection probability 
% - allowSelf - whether to allow self-connections
% - undirected - whether to make the matrix directed or not
% - ensureConnected - only return a connected matrix.
%     If this is set to true, then p must be >= 2/N(N-1) for undirected, 1/N(N-1) for directed, for the matrix to be connected.
%     Note: if allowSelf==true, then (N-1) -> N here.
% - ensureALinkForEachNode - whether to make sure that each node connects to at least one other.
%     The supplied value here is ignored if ensureConnected is set to true and an undirected graph is requested.
%     The supplied value is ignored if a directed graph is requested (because I haven't thought about whether 
%      we would want it to mean all have connections in, out or bidirectional).
%     Setting this to true will give a higher probability that the final graph is connected.
%     I don't think we can grow the graph by adding each node one at a time, giving it a link to the existing cluster,
%      then once all nodes are in adding enough extra links to get p right - I think the earlier nodes will have
%      more links - this paper describes a similar process that gave different structures to truly random graphs:
%      http://math.uchicago.edu/~shmuel/Network-course-readings/CHKNS.pdf
%
% Outputs
% - A - connectivity matrix (can be directed; A(i,j) means a link exists from i->j)
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

if (allowSelf)
    numPossibleDirectedConnections = N * N;
    withSelfSuffix = '';
else
    numPossibleDirectedConnections = N * (N - 1);
    withSelfSuffix = 'out';
end

if (nargin < 6)
    % This will be overridden if ensureConnected == true
    ensureALinkForEachNode = false;
end

if (ensureConnected)
    % If we want a connected graph, then there must be a link for each node.
    % This doesn't ensure it will be connected, but if it's not the case then the graph definitely won't be connected.
    ensureALinkForEachNode = true;
    % Now check that the user has requested enough links to make the sure the graph is connected.
    %  (Need at least N-1 links to make it undirectionally connected)
    if (undirected && (p < 2.*(N-1)./numPossibleDirectedConnections))
        error('Require p >= 2*%d/%d to have a connected undirected graph of %d nodes with%s self-connections (p=%.4f, 2*%d/%d = %.4f)\n', ...
            (N-1), numPossibleDirectedConnections, N, withSelfSuffix, p, (N-1), ...
            numPossibleDirectedConnections, 2.*(N-1)./numPossibleDirectedConnections);
    end
    if (not(undirected) && (p < 1.*(N-1)./numPossibleDirectedConnections))
        error('Require p >= 1*%d/%d to have an (undirectionally) connected directed graph of %d nodes with%s self-connections (p=%.4f, 1*%d/%d = %.4f)\n', ...
            (N-1), numPossibleDirectedConnections, N, withSelfSuffix, p, (N-1), ...
            numPossibleDirectedConnections, (N-1)./numPossibleDirectedConnections);
    end
end

while (true)
    % Generate a new network structure with connectivity depth p
    A = zeros(N);

    if (ensureALinkForEachNode)
        if (not(undirected))
            % The supplied value is ignored if a directed graph is requested (because I haven't thought about whether 
            %  we would want it to mean all have connections in, out or bidirectional).
            fprintf('Argument ensureALinkForEachNode=true is currently ignored for directed graphs!\n');
            pBetweenAll = p;
        else
            % Introduce one link for each node. In fact, this adds one more link that is strictly necessary (N-1 are necessary here)
            %  but we will reduce p for the rest of the network to take account of this. 
            % It would be foolish to call this function with the absolute minimum p anyway since it will be very difficult
            %  to randomly generate a connected network with that p.
            if (not(allowSelf) || ensureConnected)
                % Work out who to connect each node to, removing one option representing self-connections
                nodesToConnectEachTo = ceil(rand(N, 1) .* (N - 1));
                % Make the undirected connections
                for n = 1 : N
                    if (nodesToConnectEachTo(n) >= n)
                        % This means we want to connect to a node with a higher index than ourselves, so 
                        %  adjust this node index upwards (since we removed one index to represent ourself).
                        nodesToConnectEachTo(n) = nodesToConnectEachTo(n) + 1;
                    end
                    A(n, nodesToConnectEachTo(n)) = 1;
                    A(nodesToConnectEachTo(n), n) = 1;
                end
            else
                % We're allowing self-connections and not ensuring the network is connected.
                % Work out who to connect each node to
                nodesToConnectEachTo = ceil(rand(N, 1) .* N);
                % Make the undirected connections
                for n = 1 : N
                    A(n, nodesToConnectEachTo(n)) = 1;
                    A(nodesToConnectEachTo(n), n) = 1;
                end
            end
            % Now adjust the connection probability in the rest of the network for the N undirected connections made here:
            pBetweenAll = p - 2.*N./numPossibleDirectedConnections;
        end
    else
        pBetweenAll = p;
    end

    % Assign directed connections at random, but allow any initial connections made above to come through without pBetweenAll applying:
    A = A | (rand(N) < pBetweenAll);
    if (not(allowSelf))
        % Remove self-connections
        A = A .* not(eye(N));
    end

    if (undirected)
        % Just take the upper triangle and reflect it
        A = triu(A);
        A = A | A';
        symmA = A;
    else
        symmA = A | A';
    end
    
    if (ensureConnected)
        % Now check that the matrix is connected
        visited = zeros(1,N);
        visited = runDfs(symmA,N,1,visited);
        if (sum(visited) ~= N)
            % The matrix is not undirectionally connected
            fprintf('Matrix is not unidirectionally connected, trying again ...\n');
            continue;
        end
    end
    % The matrix is fine to return:
    % Debug print:
    actualP = sum(sum(A)) ./ numPossibleDirectedConnections;
    fprintf('Generated a graph for p=%.4f with p=%.4f\n', p, actualP);

    return;
end
end

