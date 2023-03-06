function A = generateNewRandomFixedDMatrix(N, p, allowSelf, undirected, ensureConnected)
%% Generate a new random matrix with a fixed indegree for each node
%
% Inputs
% - N - network size
% - p - connection probability 
% - allowSelf - whether to allow self-connections
% - undirected - whether to make the matrix directed or not
% - ensureConnected - only return a connected matrix.
%     If this is set to true, then p must be >= 2/N(N-1) for undirected, 1/N(N-1) for directed, for the matrix to be connected.
%     Note: if allowSelf==true, then (N-1) -> N here.
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
    numInLinks = floor(p * N);
else
    numPossibleDirectedConnections = N * (N - 1);
    withSelfSuffix = 'out';
    numInLinks = floor(p * (N - 1));
end

if (ensureConnected)
    % Now check that the user has requested enough links to make the sure the graph is connected.
    %  (Need at least N-1 links to make it undirectionally connected)
    if (undirected && (p < 2.*(N-1)./numPossibleDirectedConnections))
        error('Require p >= 2*%d/%d to have a connected undirected graph of %d nodes with %s self-connections (p=%.4f, 2*%d/%d = %.4f)\n', ...
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

    if (undirected)
        % For undirected network, we need to:
        % For each node, assign the remaining number of connections needed among the nodes remaining available for connection, then take that node out of the pool.
        remainingNodes = [1:N];
        connections = zeros(1,N);
        failedToAssignLinksOk = false;
        for node = 1 : N
            if (length(remainingNodes) == 0)
                break;
            end
            if (remainingNodes(1) == node)
                % We need to assign some links for this node
                if (not(allowSelf))
                    remainingNodes(1) = []; % First remove it from the remaining pile
                end
                numNodesToSelect = numInLinks - connections(node);
                if (length(remainingNodes) < numNodesToSelect)
                    % It's actually possibly that we connected up the network to this point in such a way that there are not enough options
                    %  left for us to connect this node to. So, we leave the loop and try again.
                    failedToAssignLinksOk = true;
                    break;
                end
                [randVals, indicesToSelectFromRemainingNodes] = sort(rand(1, length(remainingNodes)));
                nodesToConnectTo = remainingNodes(indicesToSelectFromRemainingNodes(1:numNodesToSelect));
                % Now make the connections:
                A(node, nodesToConnectTo) = 1;
                A(nodesToConnectTo, node) = 1;
                % And add 1 to the connections for each node that we've linked to here:
                % (for self-links, this adds an extra 1 if we were to add numNodesToSelect, but we're just going to hard 
                % code its number of links to numInLinks
                connections(nodesToConnectTo) = connections(nodesToConnectTo) + 1;
                connections(node) = numInLinks;
                % And remove them from the available pile if their work here is done ...
                indicesOfNodesAtMaxLinksInNodesWeConnectedTo = find(not(connections(nodesToConnectTo) - numInLinks));
                if (length(indicesOfNodesAtMaxLinksInNodesWeConnectedTo) > 0)
                    % Remove these nodes from the list of remaining nodes.
                    % This will also remove the current destination node, if it was left here to allow self-links.
                    remainingNodes(indicesToSelectFromRemainingNodes(indicesOfNodesAtMaxLinksInNodesWeConnectedTo)) = [];
                end
                % And make sure we've removed ourself:
                if (allowSelf && (remainingNodes(1) == node))
                    remainingNodes(1) = [];
                end
            end
        end
        if (failedToAssignLinksOk)
            % Try again
            continue;
        end
        % Sanity check
        if (allowSelf)
            % total number of links is sum of adjacency matrix + add self-links in again)
            if (N * numInLinks ~= sum(sum(A)) + diag(A))
                error('We do not have the correct number of links for each node!');
            end
        else
            % total number of links is the sum of adjacency matrix
            if (N * numInLinks ~= sum(sum(A)))
                error('We do not have the correct number of links for each node!');
            end
        end
        symmA = A;
    else
        % For directed network, we just need to randomly assign numInLinks incoming links
        if (allowSelf)
            numCandidates = N;
        else
            numCandidates = N - 1;
        end
        for node = 1 : N
            [randVals, indicesToSelect] = sort(rand(1, numCandidates));
            if (allowSelf)
                A(indicesToSelect(1:numInLinks), node) = 1;
            else
                % Adjust the source indices to account for the fact that we can't connect to ourselves
                %  so our self-index is missing here (some indices need to be bumped up)
                sources = indicesToSelect(1:numInLinks) + (indicesToSelect(1:numInLinks) >= node);
                A(sources, node) = 1;
            end
        end
        symmA = A | A';
    end

    if (ensureConnected)
        % Now check that the matrix is connected
        visited = zeros(1,N);
        visited = runDfs(symmA,N,1,visited);
        if (sum(visited) ~= N)
            % The matrix is not *undirectionally* connected
            fprintf('Matrix is not unidirectionally connected, trying again ...\n');
            continue;
        end
    end
    % The matrix is fine to return:
    % Debug print:
    actualP = sum(sum(A)) ./ numPossibleDirectedConnections;
    fprintf('Generated a graph with constant in-degree for p=%.4f with p=%.4f\n', p, actualP);

    return;
end
end

