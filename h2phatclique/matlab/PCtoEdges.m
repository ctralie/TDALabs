function [ edgelist ] = PCtoEdges( PC )
%PCTOEDGES Takes a point cloud and outputs its complete edge and vertex
%list in a format suitable for input to Cliquer (a.k.a. findcliques)
%   Detailed explanation goes here

numvert=size(PC,1);

% number of edges is n(n-1)/2
edgelist=zeros(numvert+numvert*(numvert-1)/2,3);

% vertices are specified as (i,i) with length 0
for i=1:numvert
    edgelist(i,1:2)=i*[1,1];
end

% insert all (i,j) vertex combinations, in order corresponding to
% row-by-row entries in an upper triangular matrix of pairwise distances
for i=1:numvert
    for j=(i+1):numvert
        edgelist((i*numvert)+j-i*(i+1)/2,1)=i;edgelist((i*numvert)+j-i*(i+1)/2,2)=j;
    end
end

% insert pairwise distances
edgelist((numvert+1):end,3)=pdist(PC);

% make vertex/edge numbering start from 0, so Cliquer doesn't die horribly
edgelist(:,1:2)=edgelist(:,1:2)-1;

end

