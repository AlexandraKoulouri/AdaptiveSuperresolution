function  [x_total,bool,x_cand] = AddExtraPoints_1D(Ind,x,h_u)

% x_total: total number of grid point 
% if bool = 1 then no more extra point added
% x_cand: candidate (possible) locations where the coefficients of the
% numerical problem are not zero.

bool = 0;
ExtraPoints = [];

if length(Ind)==1
     
    %there is only one non-zero coefficient
     x_extra = Set2ExtraPoints(x,Ind,h_u) ;
     ExtraPoints = [ExtraPoints x_extra]; %#ok<*AGROW>
    
else

       %Set two new grid points between a grid point with one non-zero coefficient   
        j = 1;    
        while j<=length(Ind) 
            x_extra = Set2ExtraPoints(x,Ind(j),h_u) ;
            ExtraPoints = [ExtraPoints x_extra]; %#ok<*AGROW>
            j = j+1;
        end
          
end

x_old = x(Ind);
x_cand =sort(unique([x_old(:); ExtraPoints(:)]));
if isempty(ExtraPoints) bool =1; end  
x_total = sort(unique([x(:); ExtraPoints(:)]));



function x_extra = Set2ExtraPoints(x,Ind,h_u)
      

x_extra = [];

if Ind ==1
       h = x(Ind(1)+1)-x(Ind(1));
       if (h>0.15*h_u)
          x_extra = x(Ind(1))+h/2;
       end   

elseif Ind == length(x)

       h = x(end)-x(end-1);
       if (h>0.15*h_u)
           x_extra = x(end)-h/2;
       end
          
else

       h1 = x(Ind)-x(Ind-1);
       h2 = x(Ind+1) - x(Ind);
       x_extra1 = [];
       x_extra2 = [];
       if (h1>0.15*h_u)
             x_extra1 = x(Ind)-h1/2;
       end
       if (h2>0.15*h_u)
            x_extra2 = x(Ind)+h2/2;
       end 
      x_extra = [x_extra1 x_extra2]; %#ok<AGROW>
end