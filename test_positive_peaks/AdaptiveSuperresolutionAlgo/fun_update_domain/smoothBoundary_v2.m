function [H,e,ind_poly_total] = smoothBoundary_v2(gx,gy,H)
 
[cluster_elem]= clustering(H,gx,gy,1:size(H,1));

max_cluster = max(cluster_elem(:,2));
ind_polyg_clusters =[];
for i = 1:max_cluster
    ind_cluster = find(cluster_elem(:,2)==i);
    e =boundedges([gx gy],H(cluster_elem(ind_cluster,1),:));
    ind_polygon_line = cluster_bnd_nodes_debug(e,1,gx,gy);
    ind_polyg_clusters = [ind_polyg_clusters;[ind_polygon_line(:,1) zeros(size(ind_polygon_line,1),1)+i]];
end
    


 e =[];
 ind_poly_total =[];
for cluster_iter = 1:max(ind_polyg_clusters(:,2)) %for each disconnected polygone line

    
     N = length(ind_polyg_clusters(ind_polyg_clusters(:,2)==cluster_iter,1));
     %Estimate the inner produce and the angle between the vectors created by
     %points p2-p1 and p3-p1
    ind_poly_cluster = ind_polyg_clusters(ind_polyg_clusters(:,2)==cluster_iter,1); %indicies of the nodes of the polygon line
     ind_poly_cluster = [ind_poly_cluster;ind_poly_cluster(1)];
     %coordinates of the polygon
    xq = gx(ind_poly_cluster);
    yq = gy(ind_poly_cluster);
    xq = [xq; xq(2)];
    yq = [yq; yq(2)];
    Bool_no_empty = any(inpolygon(gx(setdiff(1:length(gx),unique(ind_polyg_clusters(:,1)))),gy(setdiff(1:length(gx),unique(ind_polyg_clusters(:,1)))),xq,yq));%check is the polygone line is inside another polygone line (the polygon line which is inside another one is not updated)



    if (Bool_no_empty==1) || max(ind_polyg_clusters(:,2))==1
    
    
    %e = [e;[tempPoly_group(1:end-1) tempPoly_group(2:end)]];
    ind_poly_total = [ind_poly_total;ind_polyg_clusters(ind_polyg_clusters(:,2)==cluster_iter,:)];
    checkagain=1;
    % figure
    % triplot(H,gx,gy)
    % hold on
    % plot(xq,yq)
    while checkagain==1
        [H,e_group,ind_poly_total,checkagain,ind_poly_cluster]= FillBoundaryElements(H,gx,gy,ind_poly_total,ind_poly_cluster);
    %plot(xq,yq,'y')
    end
     
    e = [e;e_group];
end



end
