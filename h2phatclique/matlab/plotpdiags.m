function [ ] = plotpdiags( filename )
%PLOTPDIAGS Plot persistence diagrams given pairs of points with weights
%   Detailed explanation goes here

pairs = dlmread(filename);
maxdim=pairs(end,1)+1;

% filter pairs to get rid of those with pers=0 (within floating point
% error, which I'm taking to be 10^-16.)
%%%%% NO LONGER REQUIRED -- WE NOW DO THIS FILTERING IN C++ %%%%%
%pairs(abs((pairs(:,5)-pairs(:,4)))<10^(-16),:)=[];

% separate by dimension and 
seppairs=cell(1,maxdim);
infpairs=cell(1,maxdim);
maxes=zeros(1,maxdim);
for i=1:maxdim
    seppairs{i} = pairs(pairs(:,1)==(i-1) & pairs(:,3)~=-1,4:5);
    maxes(i)=max(seppairs{i}(:));
    infpairs{i} = pairs(pairs(:,1)==(i-1) & pairs(:,3)==-1,4:5);
    infpairs{i}(:,2)=1.2*max(seppairs{i}(:,2));
    seppairs{i}=[seppairs{i};infpairs{i}];
end

% for displaying on screen
for i=1:maxdim
    if size(seppairs{i},1)>0 
        figure;

        % plot birth/death pairs
        scatter(seppairs{i}(:,1),seppairs{i}(:,2),'filled');
    
        % add line x=y if dim>0 or line x=0 if dim=0
        % and a thick black line, above which we plotted inf pers pts
        hold on
        x=linspace(0,1.25*maxes(i));
        y=zeros(size(x))+1.1*maxes(i);
        if i>1
            plot(x,x);
            plot(x,y,'k','LineWidth',2);
        else
            plot(zeros(length(x)),x);
            plot([x -x],[y y],'k','LineWidth',2);
        hold off
        end
        
    axis tight    
    end
end

% for printing to files
%{
for i=1:maxdim
    hFig=figure();
        if size(seppairs{i},1)>0 
        figure;

        % plot birth/death pairs
        scatter(seppairs{i}(:,1),seppairs{i}(:,2),'filled');
    
        % add line x=y if dim>0 or line x=0 if dim=0
        % and a thick black line, above which we plotted inf pers pts
        hold on
        x=linspace(0,1.25*maxes(i));
        y=zeros(size(x))+1.1*maxes(i);
        if i>1
            plot(x,x);
            plot(x,y,'k','LineWidth',2);
        else
            plot(zeros(length(x)),x);
            plot([x -x],[y y],'k','LineWidth',2);
        hold off
        end
        
    axis tight    
    end
    print(hFig,'-dpng',fname,'-r300');
end
%}

end

