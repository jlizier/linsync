function nowVisited = runDfs(A,N,node, visited)
%
% Run a depth first search on the network A from the given node.
% This is performed recursively, and works ok on my machine up to size N = 250 -
%  I suspect it will need to be rewritten iteratively for larger network sizes.
%
% Inputs
% - A - connectivity matrix (can be directed; A(i,j) ~= 0 means a link exists from i->j)
% - N - network size
% - node - current node index to visit in [1..N]
% - visited - boolean array, storing which nodes we've visited so far.
%    For the initial external call, either don't pass this in or set to
%    zeros(1,N)
%
% Outputs
% - nowVisited - updated array of which nodes have been visited
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

if (nargin < 4)
    visited = zeros(1,N);
end

visited(node) = true;
% printf("Visiting node %d\n", node);
for i = 1 : N
    % Check if we can visit i from node
    if (A(node, i) && not(visited(i)))
        % visit i
        visited = runDfs(A,N,i, visited);
    end
end
% Return the boolean array for which nodes we've visited.
nowVisited = visited;

end
