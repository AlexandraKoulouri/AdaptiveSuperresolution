
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