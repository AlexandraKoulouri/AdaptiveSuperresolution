function [Ind_nz_cell,c_hat_cell]= fun_arrange_nzentries_subregions(NumSubRegions, SizeOfSubRegion,NonZeroInd,c_hat)

%inputs 
% total number of non zero coefficients: NonZeroInd
% estimated coefficients: c_hat
% Number of sub-regions: NumSubRegions
% Number of grid points in each sub-region:  SizeOfSubRegion

% this function returns the indices of the grid points with non-zero
% entries at each sub-region and the values of the correspoinding points

Ind_nz_cell = cell(NumSubRegions,1);  
c_hat_cell  = cell(NumSubRegions,1);

it = 1;

upper = SizeOfSubRegion{1}(1,1); 
lower = 0;


 while it<= NumSubRegions
     
     
      i = find(NonZeroInd<=upper & NonZeroInd>lower);
      c_hat_cell{it} = c_hat(NonZeroInd(i));
      Ind_nz_cell{it} = NonZeroInd(i)-lower; 
             
      if it == NumSubRegions; break; end;
      
      it = it+1;
      lower = upper;
      upper = upper + SizeOfSubRegion{it}(1,1); 
      i=[];  
     
 end
     
     