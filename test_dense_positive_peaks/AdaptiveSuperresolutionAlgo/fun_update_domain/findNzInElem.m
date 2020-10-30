function [ind_el_nz,ind_nz_nodes]=findNzInElem(El,x,y,x_nz,y_nz)

Num_nz_nodes = size(x_nz,1);
ind_nz_nodes = zeros(Num_nz_nodes,1); 
ind_el_nz = [];   

%find elements where the non-zero entries belong to
for i = 1 : Num_nz_nodes

        ind_nz_nodes(i) = find(x == x_nz(i,1) & y == y_nz(i,1),1); 
        indEltemp = find( El(:,1)==ind_nz_nodes(i) | El(:,2)==ind_nz_nodes(i) | El(:,3)==ind_nz_nodes(i));
        ind_el_nz = [ind_el_nz; indEltemp]; %#ok<AGROW>
end
ind_el_nz = unique(ind_el_nz); 