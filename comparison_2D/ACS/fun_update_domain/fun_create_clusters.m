function [El,x,y,ind_nodes_clustered,ind_inner_nodes_clustered,clustered_elem] = fun_create_clusters(x_nz,y_nz, x,y, El)
%Inputs:
%(x_nz,y_nz) nodes of the (x,y) mesh which we want to cluster 
% El (Elements): triangulation

%Outputs: 
%1. the updated elements H 
%2. the updated (x,y) nodes of size n
%3. ind_nodes_clustered [n x 2]-> first column gives the indices of the node and the second
%column in which cluster they belong to
%4. ind_inner_nodes_clustered [n1 x 2]: only the indices of the nodes which are inside the
%polygon lines and not at the edges appear in the first column and in
%second column we have in which cluster (or polygon) they belong to
%5. clustered_elem (includes the elements where at least one node
%corresponded to non-zero entry and the respect cluster that the elements
%belong to)

Num_nz_nodes = length(x_nz);
ind_nz_nodes = zeros(Num_nz_nodes,1); 
ind_el_nz = [];   

%find elements where the non-zero entries belong to
for i = 1 : Num_nz_nodes

        ind_nz_nodes(i) = find(x == x_nz(i,1) & y == y_nz(i,1),1); 
        indHtemp = find( El(:,1)==ind_nz_nodes(i) | El(:,2)==ind_nz_nodes(i) | El(:,3)==ind_nz_nodes(i));
        ind_el_nz = [ind_el_nz; indHtemp]; %#ok<AGROW>
end
ind_el_nz = unique(ind_el_nz); 
maxH = size(El,1);

%expand mesh if the non-zero entries are located at boundary nodes
[x_nz,y_nz,x,y,El]  = fun_expand_cluster(x_nz,y_nz,x,y,El,ind_el_nz,ind_nz_nodes); %#ok<ASGLU>

if maxH<size(El,1)
   ind_el_nz = [ind_el_nz;[maxH+1:size(El,1)]']; %if extra elements are added then update the indicies of the elements with non-zero entries
 
end
    
[clustered_elem, ind_nodes_clustered, ind_inner_nodes_clustered]= clustering(El,x,y,ind_el_nz); 



