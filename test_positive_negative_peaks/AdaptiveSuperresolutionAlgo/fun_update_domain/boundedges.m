function [e, Ind_BndEl,Ind_NonBndEl] =boundedges(p,t)
%BOUNDEDGES Find boundary edges from triangular mesh
%also the Indices of the elements where the boundary edges belong
%and the indices of the elements which they do not have any boundary edge
%   E=BOUNDEDGES(P,T)



ElNum = size(t,1);

edges=[t(:,[1,2]);
       t(:,[1,3]);
       t(:,[2,3])];
node3=[t(:,3);t(:,2);t(:,1)];
edges=sort(edges,2);
[foo,ix,jx]=unique(edges,'rows');

vec=histc(jx,1:max(jx));
qx=find(vec==1);
e=edges(ix(qx),:);
node3=node3(ix(qx));

Ind_BndEl = mod(ix(qx),ElNum);
ind = find(Ind_BndEl == 0);
Ind_BndEl(ind) = ElNum;
Ind_BndEl = unique(Ind_BndEl);

%xx = t(Ind_BndEl,:);
% triplot(xx,p(:,1),p(:,2),'r')

Ind_NonBndEl = setdiff(1:ElNum,Ind_BndEl);


% Orientation
v1=p(e(:,2),:)-p(e(:,1),:);
v2=p(node3,:)-p(e(:,1),:);
ix=find(v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1)>0);
e(ix,[1,2])=e(ix,[2,1]);
