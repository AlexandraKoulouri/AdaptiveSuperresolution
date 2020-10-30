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


%find elements where the non-zero entries belong to
[ind_el_nz,ind_nz_nodes]=findNzInElem(El,x,y,x_nz,y_nz);

%El=unique(El,'rows');
% figure
% triplot(El,x,y)
% hold on
% plot(x_nz,y_nz,'x')
%expand mesh if the non-zero entries are located at boundary nodes
[x_nz,y_nz,x,y,El]  = fun_expand_cluster(x_nz,y_nz,x,y,El,ind_el_nz,ind_nz_nodes); %#ok<ASGLU>

% if maxEl<size(El,1)
%    ind_el_nz = [ind_el_nz;[maxEl+1:size(El,1)]']; %if extra elements are added then update the indicies of the elements with non-zero entries
%  
% end
%find elements where the non-zero entries belong to (possibly El structure
%has been updated!)
[ind_el_nz,ind_nz_nodes]=findNzInElem(El,x,y,x_nz,y_nz);

%%triplot(El,x,y)
[clustered_elem, ind_nodes_clustered, ind_inner_nodes_clustered,El]= clustering(El,x,y,ind_el_nz); 



