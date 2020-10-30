function cm = constructColormapGray

% map_temp=colormap(hot);
% if size(map_temp,1) == 64 
% map_temp=map_temp(64:-1:1,:);
% levels = linspace(5,50,n-3);
% 
% colMap = zeros(n,3);
% %colMap(1,:) = map_temp(end,:);
% colMap(3:end-1,:) = map_temp(round(levels),:);
% colMap(end,:)     = map_temp(end-10,:);
% else
% colMap = map_temp;
cm = zeros(63,3);
cm(1:62,:) = [ones(62,1) linspace(0,1,62)'  linspace(0,1,62)'];
cm(60:64,:)   = ones(5,3);
cm(1:4,:) = [ones(4,1) zeros(4,1) zeros(4,1)];

