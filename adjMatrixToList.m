function List = adjMatrixToList(A, listForLinksIn)
%% Convert the adjacency matrix A to an adjacency list of who the connected nodes are, 
% for connections into each node iff listForLinksIn == true, else
% for connections out of each node iff listForLinksIn == false.
% Node pair is considered connected iff A(i,j) ~= 0
% 
% Input:
% - A - 2D matrix, A(i,j) means a link i->j
% - listForLinksIn (optional) - boolean, whether to list connections into (true) or
% out from each node
% 
% Output:
% - List - cell array, where A(n) is a vector of the connected node ids for n
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3


if (nargin < 2)
    listForLinksIn = false;
else

List = cell(length(A),1);

for n=1:length(A)
    if (listForLinksIn)
        % Pull out the incoming connections for n
        List{n} = find(A(:,n) > 0)';
    else
        % Pull out the outgoing connections for n
        List{n} = find(A(n,:) > 0);
    end
end
end
