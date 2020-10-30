%return the number of clusters and the boundary edges

function  [cluster, ind_nodes_clusters, ind_inner_nodes_clusters]= clustering(El,x,y,ind_nz_el)

%This function find the elements where the non-zero entries exist and
%clusters them

%Input
%triangular mesh H
%(x,y) coordinate of mesh H
%ind_nz_H: the indicies of the elements where non-zero coefficients appear

%Output
%cluster: consists of two columns, 1st column shows that index of
%element(which has at least one non-zero coefficient) and 2nd column has
%the number of the cluster that this element belongs to

%ind_nodes_custers: returns the indicies of the nodes and in which cluster
%they belong

%ind_inner_nodes_clusters: returns the indicies of the nodes (that do not
%belong to the boundary and the number of cluster that they belong to.


%this code was created by A. Koulouri, 25.4.2016

cr = [];
R =[];

ind_nz_el = sort(unique(ind_nz_el)); %indices of elements where non-zeros belongs
nodes_ind = El(ind_nz_el,:);             %indices of nodes of each element
i=1;
while i<=size(nodes_ind ,1) %draw a circle around each element and check intersections between circles
           [Rt,crt]= Circle([ x( nodes_ind(i,:)) y( nodes_ind(i,:))]); %#ok<AGROW>
           if isempty(Rt)
              El(ind_nz_el(i),:)=[];
              nodes_ind(i,:) =[];
              ind_nz_el(i) = [];
              ind_nz_el(i:end)=ind_nz_el(i:end)-1;
             % i = i-1;
           else
               i = i+1;   
               cr = [cr;crt];%circle centre 
               R = [R;Rt]; %radius for each element indHt
           end  
end

%centre of mass of the elements
 c_x = mean(x(nodes_ind)');
 c_y = mean(y(nodes_ind)');

St = [2:size(ind_nz_el,1)]; %all the possible elements 
%iter_st = 1;

cluster = zeros(size(ind_nz_el,1),2);
clusterNo = 1; %first cluster
cluster(1,1) = ind_nz_el(1);%element index H(indHt(1),:);%first elements
cluster(1,2) = 1; %cluster 1
cur_elem_ind =1;
iter = 2;

while~isempty(St)
    
   % j = 1;
       
        while (~isempty(St))
             i = length(cur_elem_ind);
              while i>=1 &&  (~isempty(St))
                      
                      Dist12 = ( c_x(cur_elem_ind(i)) - c_x(St)).^2+(c_y(cur_elem_ind(i))-c_y(St)).^2;
                      j = find(Dist12 <=1.3*min(Dist12));
                       for jiter = 1:length(j)
                        
                               if (circcirc_intersection(cr(cur_elem_ind(i),1),cr(cur_elem_ind(i),2),R(cur_elem_ind(i)),cr(St(j(jiter)),1),cr(St(j(jiter)),2),R(St(j(jiter)))))
                                %check if the closest element to S[j] and S[j] element intersect!     
 
                                    IndNode_Stj = El(ind_nz_el(St(j(jiter))),:);
                                    IndNode_ElemInd = El(ind_nz_el(cur_elem_ind(i)),:);
                                    [~,~,ff] = unique([IndNode_Stj IndNode_ElemInd]);
                                     if length(find(hist(ff)==2))==2

                                        cluster(St(j(jiter)),2) = cluster(cur_elem_ind(i),2);
                                        cluster(St(j(jiter)),1) = ind_nz_el(St(j(jiter)));
                                        %indcluster =  find( cluster(:,2) == cluster(cur_elem_ind(i),2));
                                        cur_elem_ind(iter) = St(j(jiter)); %#ok<*AGROW>
                                        iter = iter+1; 
                                        St(j(jiter)) = [];
                                        i=length(cur_elem_ind)+1;
                                        break;                                      
                                     end

                               end
                      
                       end
                     i = i-1; %go to the next element
                     
              end
                                   
            % j = j+1;
       %if you have check all the ElemInd and there is not intersection
       %create a new group
          
       if ((i==0) && ~isempty(St))
           
            clusterNo = clusterNo+1;  
            cluster(St(1),2) = clusterNo;
            cluster(St(1),1) = ind_nz_el(St(1));
            cur_elem_ind(iter) = St(1);
            iter = iter+1;
            St(1)=[];
            i = 1;
         
       end
   end
      
        
    
    
end

nodes_ind = unique(nodes_ind);
ind_nodes_clusters = zeros(length(nodes_ind),2);
ind_inner_nodes_clusters= zeros(length(nodes_ind),2);

step1 = 0;
step2 = 0;
for i = 1:max(cluster(:,2))
    ind_elem_cluster =find(cluster(:,2) == i); 
    ind_nodes_cluster = El(cluster(ind_elem_cluster,1),:);
    ind_nodes_cluster = unique(ind_nodes_cluster(:));
   
    ind_nodes_clusters( step1+1:step1+size(ind_nodes_cluster,1),1) = ind_nodes_cluster;
    ind_nodes_clusters( step1+1:step1+size(ind_nodes_cluster,1),2) = zeros(size(ind_nodes_cluster,1),1)+i;
    step1 = step1+size(ind_nodes_cluster,1);
    T = El(cluster(ind_elem_cluster,1),:);
    ap = size(x,1); 
    A = min(sparse(T(:,1),T(:,2),1,ap,ap)+sparse(T(:,2),T(:,3),1,ap,ap)+sparse(T(:,3),T(:,1),1,ap,ap),1);
    A = min(A+A',1);
    
% % this finds the boundary points, whose indexes are stored in Ibord
    B = A^2.*A==1; 
    Ibord = find(sum(B,2)>0);
    ind_inner =setdiff(ind_nodes_cluster',Ibord)'; 
    ind_inner_nodes_clusters(step2+1:step2+length(ind_inner),1) = ind_inner;
    ind_inner_nodes_clusters(step2+1:step2+length(ind_inner),2) = zeros(length(ind_inner),1)+i;
    step2 = step2+size(ind_inner,1);
   
    
end
ind_inner_nodes_clusters(ind_inner_nodes_clusters(:,2)==0,:)=[];

function  [xout,yout]= IntersectionTwoCircles(cr1,R1,cr2,R2)

 [xout,yout] = circcirc(cr1(1),cr1(2),R1,cr2(1),cr2(2),R2);
 
 %if this algorithm returns NAN then they do not interest
 %otherwise it should return the nodes of the common edge!
 
 function [R,sc] = Circle (g)
%insert 3 points i.e g  is a 2by3 matrix
p1.x = g(1,1); p2.x=g(2,1); p3.x=g(3,1);
p1.y = g(1,2); p2.y=g(2,2); p3.y=g(3,2);
       
    if (~IsPerpendicular(p1, p2, p3))                 
        [R,sc] = CalcCircle(p1, p2, p3);   
    elseif (~IsPerpendicular(p1, p3, p2))       
         [R,sc] = CalcCircle(p1, p3, p2);   
    elseif(~IsPerpendicular(p2, p1, p3))       
        [R,sc] = CalcCircle(p2, p1, p3);   
    elseif (~IsPerpendicular(p2, p3, p1))       
        [R,sc] = CalcCircle(p2, p3, p1);   
    elseif (~IsPerpendicular(p3, p2, p1))       
        [R,sc] = CalcCircle(p3, p2, p1);   
     elseif (~IsPerpendicular(p3, p1, p2))       
        [R,sc] = CalcCircle(p3, p1, p2);   
    else  
        %TRACE("\nThe three pts are perpendicular to axis\n");
       R = []; sc = [];
   end

function Bool = IsPerpendicular(p1, p2, p3)
% Check the given points are perpendicular to x or y axis
eps = 0.00000001;
yDelta_a= abs(p2.y - p1.y);
xDelta_a= abs(p2.x - p1.x);
yDelta_b= abs(p3.y - p2.y);
xDelta_b= abs(p3.x - p2.x);

if (xDelta_a <= eps && yDelta_b <= eps)
    %The points are pependicular and parallel to x-y axis
    Bool = 0;
elseif(yDelta_a <= eps || yDelta_b <= eps || xDelta_a <= eps || xDelta_b <= eps)
    Bool = 1;
else
    Bool = 0;
end
   
function [R, sc] = CalcCircle( p1, p2, p3)
%Function CalcCircle: Estimate the circle: Radius and centre given 3 points.

yDelta_a= p2.y - p1.y;
xDelta_a= p2.x - p1.x;
yDelta_b= p3.y - p2.y;
xDelta_b= p3.x - p2.x;
   
 if(abs(xDelta_a) <= eps && abs(yDelta_b) <= eps)
     sc(1,1)= 0.5*(p2.x + p3.x);
     sc(1,2)= 0.5*(p1.y + p2.y);
       
 else
     aSlope=yDelta_a/xDelta_a;
    bSlope=yDelta_b/xDelta_b;
    % calculate center
    sc(1,1) = (aSlope*bSlope*(p1.y - p3.y) + bSlope*(p1.x + p2.x)- aSlope*(p2.x+p3.x) )/(2* (bSlope-aSlope) );
    sc(1,2) = -1*(sc(1,1) - (p1.x+p2.x)/2)/aSlope +  (p1.y+p2.y)/2;
 end
 R =  sqrt((sc(1,1)-p1.x) * (sc(1,1)-p1.x) + (sc(1,2)-p1.y) * (sc(1,2)-p1.y)); %Radius...
% 
% function d = distanceP(x,y)
%         d = sqrt((x(1)-x(2)).^2+(y(1)-y(2)).^2);
%  % [e, Ind_BndEl,Ind_NonBndEl] =boundedges(p,t)
% 
% function Bool = InsideTringle(x,y,xtri,ytri)
% b = [ x-xtri(3);y-ytri(3)];
% A = [xtri(1)-xtri(3) xtri(2) - xtri(3);ytri(1)-ytri(3) ytri(2)-ytri(3)];
% Bool =0;    
% l = A\b;
% if l(1)>=0 && l(2)>=0 && l(1)+l(2)<=1
%     Bool=1;
% end
