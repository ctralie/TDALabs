function [ ] = makedimacs( edgelist, numvertices, outfile )
%MAKEDIMACS This function takes a weighted edge list (e.g. outputted from
%ctree), and creates a file that can be read using phatclique.
%   Detailed explanation goes here

fileID=fopen(outfile,'w');

edges=edgelist(edgelist(:,3)~=0,:);
numedges=size(edges,1);
fprintf(fileID,'p edge %i %i\n',numvertices,numedges);

for i=1:numedges
    fprintf(fileID,'e %i %i %.16f\n',edges(i,1),edges(i,2),edges(i,3));
end

fclose(fileID);

end

