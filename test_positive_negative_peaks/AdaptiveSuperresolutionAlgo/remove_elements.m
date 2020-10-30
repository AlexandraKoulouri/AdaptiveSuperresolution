function h = remove_elements(x,y,h)
%Refine mesh

s1 = size(h,1);

for i = 1:s1
    xtri = x(h(i,:));
    ytri = y(h(i,:));
     s = side_length(xtri,ytri);
     p(i) = 1/2*sum(s); %triangle surface
     A2(i) =sqrt(p(i)*(p(i)-s(1))*(p(i)-s(2))*(p(i)-s(3))); 

end
A2n = A2/max(A2);

%remove elements with very long edges and small surface
%remInd  = [];
remInd = [find(p>median(p)*2)];
% figure
% triplot(h,x,y);
% hold on 
% for i=1:length(remInd)
%     hold on
% patch(x(h(remInd(i),:)),y(h(remInd(i),:)),[227 227 227]./255,'EdgeColor',[160 160 160]./255)
% end

h(remInd,:) = [];
A2n(remInd) = [];
%Estimate the elements with the highest and smallest surface and remove
%them
remInd = find(A2n>median(A2n)*10 | A2n<median(A2n)*0.1);
h(remInd,:) = [];

% %check for possible holes in the mesh
% e = boundedges([x y],h);
% Ind_Polyg_Group = GroupEdgePoints_v2(e,1);
% 
% if max(Ind_Polyg_Group(:,2))>=2 %if this is true then there are holes in the mesh
%    %find longest polygon first
%    PolLen = zeros(max(Ind_Polyg_Group(:,2)),1);
%    for i = 1:max(Ind_Polyg_Group(:,2))
%        PolLen(i) = Poly_Length(Ind_Polyg_Group(Ind_Polyg_Group(:,2)==i,1),x,y);
%    end
%    IndMaxPolLen = find(PolLen == max(PolLen));
%    
%    HoleInd = setdiff(1:max(Ind_Polyg_Group(:,2)),IndMaxPolLen);
%    for i =1:length(HoleInd)
%        PolyInd = Ind_Polyg_Group(Ind_Polyg_Group(:,2)==HoleInd(i),1);
%        %cc = mean([xx(PolyInd),yy(PolyInd)]);
%        H1 = delaunay(x(PolyInd),y(PolyInd));
%        h = [h;PolyInd(H1)];
%    end
%        


%end
%combine the two meshes

% triplot(h,x,y);
% hold on
% %draw remaining elements
% for i=1:length(h)
%     hold on
% patch(x(h(i,:)),y(h(i,:)),[100 100 100]./255,'EdgeColor',[160 160 160]./255)
% end

function PolLeng = Poly_Length(e,x,y);
PolLeng = 0;
for i = 1:length(e)-1
    PolLeng = PolLeng+(sqrt((x(e(i))-x(e(i+1)).^2+(y(e(i))-y(e(i+1)).^2))));
end



function s = side_length(x,y)
s = zeros(3,1);
s(1) = sqrt((x(1)-x(2)).^2+(y(1)-y(2)).^2);
s(2) = sqrt((x(2)-x(3)).^2+(y(2)-y(3)).^2);
s(3) = sqrt((x(1)-x(3)).^2+(y(1)-y(3)).^2);