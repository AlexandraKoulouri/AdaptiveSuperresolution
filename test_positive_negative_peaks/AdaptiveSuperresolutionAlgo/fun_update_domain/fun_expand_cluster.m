function  [x_nz,y_nz,x,y,El,ed,ind_nz_el,ind_nz_nodes]  = fun_expand_cluster(x_nz,y_nz,x,y,El,ind_nz_el,ind_nz_nodes) 
%this function updates/expands the triangulation if (x_nz,y_nz) are  boundary points by
%adding extra nodes to the existing triangular mesh


%check first the case where there are only three nodes and they are colinear!!!
if size(El,1) == 1 %% && Num_P==3
             
    [h,s,Ar] = Triangle(x,y); 
    indt=[]; %#ok<NASGU>
    
      if Ar<=eps
          
                node_1 = [x(1) y(1)];
                node_2 = [x(2) y(2)];
                node_3 = [x(3) y(3)];
                inds = find(s == max(s));
                rn = (node_2-node_3)/s(2);
                rnp= [-rn(2), rn(1)];
            if (inds == 1) %p3 is in the middle

                Pt =(node_3+rnp/2*(s(3)+s(2))); %(p3-rnp/2*(s(3)+s(2)))'];
                indt = 3;

           elseif (inds == 2)%p1 is in the middle
                 Pt =(node_1+rnp/2*(s(3)+s(1)));% (p1-rnp/2*(s(3)+s(1)))'];
                 indt = 1;
           else %p2 is in the middlex
                 Pt =(node_2+rnp/2*(s(1)+s(2))) ;
                 indt = 2;
%               
                
            end
            for i = 1 : length(x_nz)
                IndP_shift = find(x(indt) == x_nz(i,1) & y(indt) == y_nz(i,1),1); 
            end
            x(indt) = Pt(1,1); %substitute one node(collinear) with a shifted one;
            y(indt) = Pt(1,2) ;
            if ~isempty(IndP_shift)
                x_nz(IndP_shift) =Pt(1,1);
                y_nz(IndP_shift) =Pt(1,2);
            end
                  
      end
      
    
end


%find the boundary edges
ed = [];
ed =boundedges([x y],El); 


%this functions corrects the mesh (if there is a hole or a missing boundary
%elements)
[El,ed,ind_ed_clusters] = smoothBoundary(x,y,ed,El);
ind_nz_nodes_on_boundary = find(ismember(ind_nz_nodes,unique(ed(:)))==1); %find the indices of the nodes corresponding to non-zero entries if there are lie on the boundary of the domain

% for i =1:length(ed)
%    plot(x(ed(i,:)),y(ed(i,:)),'y-') 
% end
 

%Sort ind_nz_nodes_on_boundary (so you go one after the other node)
cx = mean(x);
cy = mean(y);
ang = mod(360/(2*pi)*atan2((y(ind_nz_nodes(ind_nz_nodes_on_boundary))-cy),(x(ind_nz_nodes(ind_nz_nodes_on_boundary))-cx)),360);

[~,indsorted] = sort(ang);
ind_nz_nodes_on_boundary = ind_nz_nodes_on_boundary(indsorted);


    
    
%Extend the domain if nodes (with non-zero entries) are on the boundary of
%the cluster
if ~isempty(ind_nz_nodes_on_boundary) 
        %find the edges which are between point
        %x(ind_nz_nodes(ind_nz_nodes_on_boundary(1)))
        indedg2 =[];
        indedg1 =[];
        cur_ind_nz_node = ind_nz_nodes(ind_nz_nodes_on_boundary(1));
        indedg1 = find(ed(:,1) == cur_ind_nz_node);
        indedg2 = find(ed(:,2)== cur_ind_nz_node);
            
      if length(indedg2)==2 && length(indedg1)==2 %remove points which connect two area..

                delInd = find(x(ind_nz_nodes(ind_nz_nodes_on_boundary(1)))==x_nz & y(ind_nz_nodes(ind_nz_nodes_on_boundary(1)))==y_nz);
                x_nz(delInd) = [];
                y_nz(delInd) = [];
                ind_nz_nodes(ind_nz_nodes_on_boundary(1))=[];
      else

                cur_ind_nodes = zeros(3,1);
                cur_ind_nodes(1) = cur_ind_nz_node;
                node_1 = [x(cur_ind_nz_node), y(cur_ind_nz_node)];
                %take the first non empty
                if ~isempty(indedg1)
                    node_2 = [x(ed(indedg1(1),2)) y(ed(indedg1(1),2)) ];
                    cur_ind_nodes(2) = ed(indedg1(1),2);
                    if length(indedg1)==2 && isempty(indedg2)
                       node_3 = [x(ed(indedg1(2),2)) y(ed(indedg1(1),2)) ];
                       cur_ind_nodes(3) = ed(indedg1(2),2);
                    elseif length(indedg1)==1 && length(indedg2)==1
                        node_3 = [x(ed(indedg2(1),1)) y(ed(indedg2(1),1)) ];
                        cur_ind_nodes(3)= ed(indedg2(1),1);
                    elseif length(indedg1)==1 && isempty(indedg2)==0
                       dist('the code has an error')
                    else
                        %this case appears when two clusters are connected
                        %with a single node (very rare case because we performe smoothing to the meshes/clusters) 
                        cur_ind_edg = [ed(indedg1(2),2); ed(indedg2(:),1)];
                        ind_cluster_node2 = find(ind_ed_clusters(:,1) == ed(indedg1(1),2));
                        thiscluster = ind_ed_clusters(ind_cluster_node2,2); % find the cluster where node_2 belongs %node_2 = [x(ed(indedg1(1),2)) y(ed(indedg1(1),2)) ];
                        
                        Ang = 360;
                        for ii = 1:length(cur_ind_edg) %the find  node_3 (give node_1 and node_2) that does not belong to the same cluster/group
                            ind_check_cluster = find(ind_ed_clusters(:,1) ==  cur_ind_edg(ii));
                            
                            if (ind_ed_clusters(ind_check_cluster,2) ~= thiscluster)
                                 
                                  node_temp = [x(cur_ind_edg(ii)) y(cur_ind_edg(ii))];
                                  [~,s,~] = Triangle([node_1(1) node_2(1) node_temp(1)],[node_1(2) node_2(2) node_temp(2)]); 
                               
                                  if Ang>AngleVec(node_1,node_2,node_temp,s)
                                      node_3 = node_temp;
                                      %Ang = AngleVec(node_1,node_2,node_3,s);
                                      cur_ind_nodes(3) = cur_ind_edg(ii);
                                      break;
                                  end
                            end
                        end 
                    end


               else
                     node_2 = [x(ed(indedg2(1),1)) y(ed(indedg2(1),1)) ];
                     cur_ind_nodes(2) = ed(indedg2(1),1); cur_ind_nodes(3) = ed(indedg2(2),1);
                     node_3 = [x(ed(indedg2(2),1)) y(ed(indedg2(2),1)) ];


                end
                node_1 = [x(cur_ind_nz_node), y(cur_ind_nz_node)];
                cur_ind_nodes = cur_ind_nodes';
   
                [x,y,El,ind_nz_el] =add_elements(node_1,node_2,node_3,El,x,y,cur_ind_nodes,ind_nz_el,ind_ed_clusters); 
            end
            %continue expansion
            [x_nz,y_nz,x,y,El,ed,ind_nz_el,ind_nz_nodes]  = fun_expand_cluster(x_nz,y_nz,x,y,El, ind_nz_el,ind_nz_nodes);

%   end
    end




function [x,y,El,Elind] = add_elements(p1,p2,p3,El,x,y,indp,Elind,Ind_Polyg_Group)
%Mesh H-(x,y) is expanded around the edges p2-p1 and p3-p1
%The extra elements depend on the shape of edges p2-p1 and p3-p1
    
    
%coordinates of the polygon
 xq = x(Ind_Polyg_Group(:,1),1);
 yq = y(Ind_Polyg_Group(:,1),1);

 %Estimate the heights, circumference and surface of the triangle p1p2p3
[h,s,Ar] = Triangle([p1(1) p2(1) p3(1)],[p1(2) p2(2) p3(2)]); 
 dd = h(2);
 
 %Estimate the inner produce and the angle between the vectors created by
 %points p2-p1 and p3-p1
  Angle = AngleVec(p1,p2,p3,s);

 %check if the edges p1-p2 and p1-p3 are convex (Bool=1) or concave (Bool=0)
 convex_true = CheckConvexity(p1, p2, p3, xq,yq);

 %find projection of p1 on the line p2p3
 ph = projectionOnLine(p1,p2,p3,s); 
 tol = 1e-3; %tolerance to check where points are located (in our out of the polygone line)
 maxH = max(El(:));
if (collinear([p1; p2; p3;],1e-7)) %(Ar<=1e-7)) || colinear points or almost collinear
     %check if they form an element and discard it
      inds = find(s == max(s));
      rn = (p2-p3)/s(2); 
      rnp= [-rn(2), rn(1)];
      rnp = rnp./norm(rnp);
      %this code places two points Pt both sides of the line (then we check
      %which of the two is in the polygone line!)
    if (inds == 1) %p3 is in the middle
             
        Pt =[(p3+rnp*tol*(s(3)+s(2)))' (p3-rnp*tol*(s(3)+s(2)))'];
         p=p3;
         El = [El;indp(1) maxH+1 indp(3);indp(2) maxH+1 indp(3)];
         Elind =[Elind;size(El,1)-1;size(El,1)] ;  
    elseif (inds == 2)%p1 is in the middle
         Pt =[(p1+rnp*tol*(s(3)+s(1)))' (p1-rnp*tol*(s(3)+s(1)))'];
         El = [El;indp(1) maxH+1  indp(2);indp(3) maxH+1  indp(1);];
         p=p1; 
         Elind =[Elind;size(El,1)-1;size(El,1)] ; 
    else %p2 is in the middle
         Pt =[(p2+rnp*tol*(s(1)+s(2)))' (p2-rnp*tol*(s(1)+s(2)))'];
          p=p2;
          
          El = [El;indp(1) maxH+1 indp(2);indp(3) maxH+1 indp(2)];
          Elind =[Elind;size(El,1)-1;size(El,1)] ;
    end
      
     %check if the point is inside the the polygone line
     indPt = find(inpolygon(Pt(1,:),Pt(2,:),xq,yq)==0);
    if isempty(indPt)
           keyboard;  %remove the point from the list This has to be added!
           
           
    end
   if indPt == 1 
         P = (p+rnp*mean(s(setdiff([1:3],inds))))';
         x = [x;P(1,:)'];
         y = [y;P(2,:)'];
     elseif indPt == 2
         P = (p-rnp*mean(s(setdiff([1:3],inds))))';
         x = [x;P(1,:)'];
         y = [y;P(2,:)'];
   
   else
      
        P = (p+rnp*mean(s(setdiff([1:3],inds))))';
        
         
   end   
     
         

     
elseif convex_true == 0 %concave options

     if  Angle<=155 && Angle>30
          El = [El;indp]; %create an element from p1 and p2 and p3
          P = [];  
          Elind =[Elind;size(El,1)] ;
     elseif Angle>155
         
         rn = (p2-p3)/s(2);
         rnp= [-rn(2), rn(1)];
         %rnp = rnp/norm(rnp);
         Pt =[(p1+rnp*tol*(s(1)+s(3)))' (p1-rnp*tol*(s(1)+s(3)))'];
         indPt = []; %#ok<NASGU>
         indPt = find(inpolygon(Pt(1,:),Pt(2,:),xq,yq)==0);
         if isempty(indPt)
             keyboard;
         end
         %if length(indPt)==1 
             if indPt == 1
                 P = (p1+rnp*0.3*(s(1)+s(3)))';
                 
             elseif indPt == 2
                 P = (p1-rnp*0.3*(s(1)+s(3)))';
                 
             else
                 PA = (p1+rnp*0.3*(s(1)+s(3)));
                 PB =  (p1-rnp*0.3*(s(1)+s(3)));
                 dA = max(distanceP([p2(1) PA(1)]',[p2(2) PA(2)]'), distanceP([p3(1) PA(1)]',[p3(2) PA(2)]') );
                 dB = max(distanceP([p2(1) PB(1)]',[p2(2) PB(2)]'), distanceP([p3(1) PB(1)]',[p3(2) PB(2)]') );
                 if (dA<dB) 
                     P = PA(:);
                 else
                     P = PB(:);
                 end

             end    
          x = [x;P(1,:)];
          y = [y;P(2,:)];
      
          El = [El;indp(1) maxH+1 indp(2);indp(3) maxH+1 indp(1)];
          Elind =[Elind;size(El,1)-1;size(El,1)] ;

     elseif Angle<= 30

         %Move p1 close to p2p3 edge to change the shape of the element
         rn = (p2-p3)/s(2);
%          rn = 0.5*(p2-p3)+p3;
%          P = 0.4*(p1-rn)+rn;
%          P = P(:);
         rnp= [-rn(2), rn(1)];
         P =(p1-rnp*0.5*dd)';
         x(indp(1)) = P(1,1);
         y(indp(1)) = P(2,1);
         P = [];
         El = [El;indp]; %add a single element which is created from p1 and p2 and p3
         Elind =[Elind;size(El,1)] ;
       end
elseif (convex_true== 1) %convex option       

         if Angle>=75 && Angle<=120
          %P = [PointAntisymmetric(p1,p2,0.3*(s(1)+s(3)))'  PointAntisymmetric(p1,p3,0.3*(s(1)+s(3)))'];
           P = [PointAntisymmetric(p1,p2,0.8*min(s(1),s(3)))'  PointAntisymmetric(p1,p3,0.8*min(s(1),s(3)))'];
         
          El = [El; 
              maxH+1 indp(3) indp(1);
              indp(1) maxH+1 maxH+2;
              maxH+2 indp(2)  indp(1)
             ];  
         Elind =[Elind;size(El,1)-2;size(El,1)-1;size(El,1)] ;
         x = [x;P(1,:)'];
         y = [y;P(2,:)'];

     elseif Angle>120
        
         rn = (p2-p3)/s(2);
         rnp= [-rn(2), rn(1)];

         Pt =[(p1+rnp*tol*(s(1)+s(3)))' (p1-rnp*tol*(s(1)+s(3)))'];
         indPt = find(inpolygon(Pt(1,:),Pt(2,:),xq,yq)==0);
         if isempty(indPt)
             keyboard;
         end
         if indPt == 1
             P = (p1+rnp*0.35*(s(1)+s(3)))';
         elseif indPt == 2
             P = (p1-rnp*0.35*(s(1)+s(3)))';

         else
             %check this spot!
              P = (p1+rnp*0.35*(s(1)+s(3)))';
             
         end   

        El = [El;indp(1) maxH+1 indp(2);indp(3) maxH+1 indp(1)];
        Elind =[Elind;size(El,1)-1;size(El,1)]; 
          x = [x;P(1,:)'];
         y = [y;P(2,:)'];

     elseif Angle<75
         %it may require improvement this here!
         P = PointAntisymmetric(p1,ph,0.7*dd)';
         rn = (P-p1')/sqrt((P(1,1)-p1(1,1)).^2+(P(2,1)-p1(1,2)).^2);
         rnp= [-rn(2), rn(1)];
         Pt =[(p1+rnp*0.7*s(2))' (p1-rnp*0.7*s(2))']; 
         indPt = find(inpolygon(Pt(1,:),Pt(2,:),xq,yq)==0); %#ok<EFIND>
         if isempty(indPt)
             keyboard;   
         end  
         %[ sqrt((Pt(1,1)-p2(1,1)).^2+(Pt(2,1)-p2(1,2)).^2),sqrt((Pt(1,1)-p2(1,1)).^2+(Pt(2,1)-p2(1,2)).^2)

         %keyboard;
         P = [P Pt];
         El = [El; 
              maxH+1  maxH+2  indp(1);
              indp(1) maxH+2  indp(3);
              maxH+1  indp(1) maxH+3;
              indp(1) maxH+3  indp(2);
                           ];  
         Elind =[Elind;size(El,1)-3;size(El,1)-2;size(El,1)-1;size(El,1)]; 
         x = [x;P(1,:)'];
         y = [y;P(2,:)'];

     end
   
end 
   
function d = distanceP(x,y)
        d = sqrt((x(1,:)-x(2,:)).^2+(y(1,:)-y(2,:)).^2);
 
 function [h,s,A] = Triangle(x,y)
s = zeros(3,1);
s(1) = sqrt((x(1)-x(2)).^2+(y(1)-y(2)).^2);
s(2) = sqrt((x(2)-x(3)).^2+(y(2)-y(3)).^2);
s(3) = sqrt((x(1)-x(3)).^2+(y(1)-y(3)).^2);
p = 1/2*sum(s); %triangle surface
A =sqrt(p*(p-s(1))*(p-s(2))*(p-s(3)));
h = 2*A./s;%triangle hights

function p = PointAntisymmetric(p1,p2,leng)
%find a point antisymmetric to an edge

ang =  mod(360+180/pi*atan2(p2(1,2)-p1(1,2),p2(1,1)-p1(1,1)),360);
%leng = sqrt((p2(1,2)-p1(1,2))^2+(p2(1,1)-p1(1,1))^2);

p = leng*[cosd(ang+180) sind(ang+180)]+p1;

function p = projectionOnLine(p1,p2,p3,s)
    %estimate a point normal to the line p2-p3
    rn = (p3'-p2')./s(2);
    p  = dotProd(rn,(p1-p2)')*rn+p2';
    
    p = p';


    
function d = dotProd(p1,p2)   
    d = p1(1)*p2(1)+p1(2)*p2(2);
 
function Bool = CheckConvexity(p1, p2, p3, xq, yq)
%(xq,yq) is the polygon line around the mesh
    
%A point p inside a triangle
l1 = 0.4;
l2 = 0.3;
p = l1*p1'+l2*p2'+(1-l1-l2)*p3';

%If Bool is 1 then it is a convex otherwise it is concave
Bool = inpolygon(p(1),p(2),xq,yq);  
    
function Angle = AngleVec(p1,p2,p3,s)
 %estimate the angle between two vectors (p1p2) and (p1p3).
 innerProd = (p2(1,2)-p1(1,2))*(p3(1,2)-p1(1,2))+(p2(1,1)-p1(1,1))*(p3(1,1)-p1(1,1));
 ratio = innerProd/(s(1)*s(3));
 Angle = 180/pi*acos(ratio); 

 

    
 function [H,e,ind_poly_total] = smoothBoundary(gx,gy,e,H)
 
 %create the polygon line around the triangular mesh
 ind_polyg_groups = cluster_bnd_nodes(e,1);
 e =[];
 ind_poly_total =[];
for groupiter = 1:max(ind_polyg_groups(:,2)) %for each disconnected polygone line

    
     N = length(ind_polyg_groups(ind_polyg_groups(:,2)==groupiter,1));
     %Estimate the inner produce and the angle between the vectors created by
     %points p2-p1 and p3-p1
    ind_poly_group = ind_polyg_groups(ind_polyg_groups(:,2)==groupiter,1); %indicies of the nodes of the polygon line
     %coordinates of the polygon
    xq = gx(ind_poly_group);
    yq = gy(ind_poly_group);
    xq = [xq; xq(2)];
    yq = [yq; yq(2)];
    Bool_no_empty = any(inpolygon(gx(setdiff(1:length(gx),unique(ind_polyg_groups(:,1)))),gy(setdiff(1:length(gx),unique(ind_polyg_groups(:,1)))),xq,yq));%check is the polygone line is inside another polygone line (the polygon line which is inside another one is not updated)



    if (Bool_no_empty==1) || max(ind_polyg_groups(:,2))==1
    
    
    %e = [e;[tempPoly_group(1:end-1) tempPoly_group(2:end)]];
    ind_poly_total = [ind_poly_total;ind_polyg_groups(ind_polyg_groups(:,2)==groupiter,:)];
    checkagain=1;
    % figure
    % triplot(H,gx,gy)
    % hold on
    % plot(xq,yq)
    while checkagain==1
        [H,e_group,ind_poly_total,checkagain,ind_poly_group]= FillBoundaryElements(H,gx,gy,ind_poly_total,ind_poly_group);
    %plot(xq,yq,'y')
    end
     
    e = [e;e_group];
end



end




function [H,e,ind_poly_total,checkagain,ind_poly]= FillBoundaryElements(H,gx,gy,ind_poly_total,ind_poly)
     
xq = gx(ind_poly);
yq = gy(ind_poly);
xq = [xq; xq(2)];
yq = [yq; yq(2)];
    
ind_poly_1 = [ind_poly;ind_poly(2,:)];     

checkagain = 0;

%check if the edges p1-p2 and p1-p3 are convex (Bool=1) or concave (Bool=0)
i=2;
while i<=length(xq)-1
     
     p1 = [xq(i) yq(i)];
     p2 = [xq(i-1) yq(i-1)];
     p3 = [xq(i+1) yq(i+1)];
     Bool = CheckConvexity(p1, p2, p3, xq,yq);

     if Bool==0 %if concave then fill the missing element
              %then re-align it!
              %Estimate the heights, circumference and surface of the triangle p1p2p3
            [h,s,Ar] = Triangle([p1(1) p2(1) p3(1)],[p1(2) p2(2) p3(2)]); 

        %      %Estimate the inner produce and the angle between the vectors created by
             %points p2-p1 and p3-p1
             Angle = AngleVec(p1,p2,p3,s);
             if Angle<90
%                      plot(p1(1,1),p1(1,2),'o')
%                 plot(p2(1,1),p2(1,2),'x')
%                 plot(p3(1,1),p3(1,2),'x')
                 indp = [ind_poly_1(i,1) ind_poly_1(i-1,1) ind_poly_1(i+1,1)];
                 H = [H;indp]; %add a single element which is created from p1 and p2 and p3
                 %remove a boundary and update with a new one!since new
                 %element was created!
                 if ~isempty(intersect(find(ind_poly==ind_poly_1(i)),1))
                     
                     ind_poly_total(ind_poly_total(:,1)== ind_poly_1(i),:)=[];
                     ind_poly_total=[ind_poly_total;ind_poly_total(1,:)];
                     ind_poly(ind_poly==ind_poly_1(i))=[];
                     ind_poly = [ind_poly;ind_poly(1)];
                     ind_poly_1 = [ind_poly;ind_poly(2,:)]; 
                     
                 else
                     
                     ind_poly_total(ind_poly_total(:,1)== ind_poly_1(i),:)=[];
                     
                     ind_poly(ind_poly==ind_poly_1(i))=[];
                     ind_poly_1 = [ind_poly;ind_poly(2,:)];  
                 end
                 xq = gx(ind_poly);
                 yq = gy(ind_poly);
                 xq = [xq; xq(2)]; %#ok<AGROW>
                 yq = [yq; yq(2)]; %#ok<AGROW>
               %  i=2;
                 %tempPoly_group1(i)=[]; %remove p1 for the polygone line!
                % e()=[]; %remove the two edges and replace them with a single edge
                 checkagain = 1; 
%                  figure
%                  i
%               plot(xq,yq,'y')
             else
              i=i+1;
             end
     else
         i= i+1;
     end
end
    %update edges!
 e = [ind_poly(1:end-1) ind_poly(2:end)];
%       hold on
% for i =1:length(e)
%    plot(gx(e(i,:)),gy(e(i,:)),'r-') 
% end
% e;
 
