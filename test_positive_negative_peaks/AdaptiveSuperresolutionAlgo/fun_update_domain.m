function [xnew, ynew, total_add_true,Hnew]= fun_update_domain(ind_nz,x,y,mindist,H)

%Add extra points in the vicinity of the nodes (x,y), which have non zero
%weights/coefficients, given by ind_nz
%In this code, the new points are the centroides of the largest elements
%only (the small elements are not updated...this allow better stability aka less grid points and avoid overfitting problems)
%defined by
%mindist: is the minimum distance between two points that two nodes can
%have
%H: is the triangulation of the 2D domain (x,y) (organized in cells which
%correspond to the different clusters) 
%The second column of H tells if a cluster has reached its minumum size
%1: not yet, 2: yes, 0: there is not element
%NumCluster: is the number of sets which the points [x(Ind) y(Ind)] can be grouped 
%the number of sets depends on the number of convex area that can be
%created based on the points of the nonzero weights denoted by [x(Ind) y(Ind)].

%this code was created by A. Koulouri July 2016
%Upadated in 11.2.2020


total_add_true = 0; %this index becomes one if at least one new node is included in the domain

true_add = 0;%#ok<NASGU> %True or false  index ->if trueAdd = 0; then no more extra points are added in the existing clusters
cluster_iter = 0;    
incr_cluster = 1;%add a new cluster index

while cluster_iter<size(H,1)
    cluster_iter = cluster_iter+1; 
    true_add = 0;
    if ~isempty(ind_nz{cluster_iter,1})  
         
        if ~isempty(H{cluster_iter,1}) && H{cluster_iter,2}(1,1)~=0
                [ H_up,x_up, y_up, ind_nodes_in_clusters,ind_nodes_inside_clusters,clusters] =  fun_create_clusters(x{cluster_iter,1}(ind_nz{cluster_iter}),y{cluster_iter,1}(ind_nz{cluster_iter}),x{cluster_iter,1}(:,1),y{cluster_iter,1}(:,1),H{cluster_iter,1});
                ind_nodes_in_clusters = unique(ind_nodes_in_clusters,'rows','stable'); 
                numofclusters = max(clusters(:,2));

                
               for clusterno = 1: numofclusters %add new nodes to each of the clusters and create new independant-disjoint meshes

                     true_add =0;
                     [El_div,x_div,y_div] = fun_div_mesh(x_up,y_up,ind_nodes_in_clusters,H_up,clusterno,clusters); %define mesh H_up based on the created cluters

                     if ~isempty(El_div) 

                       [extra_nodes,true_add] = fun_add_nodes(x_div,y_div,El_div,mindist);
                      
                       ind_loc = unique(ind_nodes_inside_clusters(ind_nodes_inside_clusters(:,2)==clusterno,1)); %add only the old nodes that they corresponded to non-zero entries
                       %ind_loc = unique(ind_nodes_in_clusters(ind_nodes_in_clusters(:,2)==clusterno,1));
                       if true_add == 1
                          total_add_true = 1;
                          nodes = unique([[ x_up(ind_loc);extra_nodes(:,1)] [ y_up(ind_loc);extra_nodes(:,2)]],'rows','stable');
                     
                          
                       else
                          
                           nodes = unique([ x_up(ind_loc,1)  y_up(ind_loc,1)],'rows','stable');

                       end
                       if size(nodes,1)>2
                           Tri = delaunayTriangulation(nodes(:,1),nodes(:,2)); 
                           Tri = Tri.ConnectivityList;
                            %ed = [];
                            ed =boundedges([x_div y_div],El_div); 
                            ind_poly = cluster_bnd_nodes(ed,1);
                            xq = x_div(ind_poly(:,1));
                            yq = y_div(ind_poly(:,1));
                            
                           [mid_nodes,El_ind] = fun_middle_vert(nodes, Tri);
                           [El_ind_rem] =fun_elem_outside_polygon([xq(:) yq(:)] ,mid_nodes,El_ind);
                           Tri(El_ind_rem,:)=[];
                           
                           %correct shape of the triangular mesh
                           if size(Tri,1)>20
                                Tri = remove_elements(nodes(:,1),nodes(:,2),Tri);% check and remove some very small and very large elements
                            end

                         %  Tri = alphaTriangulation(alphaShape(nodes(:,1),nodes(:,2),1);%,'HoleThreshold ',AvgDist);

                           Hnew{incr_cluster,1} =  Tri;%#ok<*AGROW> %.Triangulation; %#ok<*AGROW>
                           Hnew{incr_cluster,2} =  true_add; %A new mesh has been created!
                            
                       else
                           
                            if size(nodes,1)==1 %create a small circle around this point!
                                px = nodes(1); py=nodes(2);
%                                 Tri=[]; nodes =[];
                                nodes =[];
                                nodes(:,1) =[px; px+1.5*mindist*cosd([0:45:366-45])']';
                                nodes(:,2) = [py; py+1.5*mindist*sind([0:45:366-45])']';
                                Hnew{incr_cluster,2} = 2; %the smallest possible mesh has been created!
                                Tri = delaunayTriangulation(nodes(:,1),nodes(:,2));
                                Hnew{incr_cluster,1} = Tri.ConnectivityList;
                                                         
                       
                          else

                            Hnew{incr_cluster,1} = [];
                            Hnew{incr_cluster,2} = 0; %no  mesh has been created in this round!
                          end
                       end

                       xnew{incr_cluster,1} = nodes(:,1); %#ok<NASGU>
                       ynew{incr_cluster,1} = nodes(:,2); %#ok<NASGU>
                       incr_cluster = incr_cluster +1;
                     
%                     
                     end
                     

                end
       
       elseif isempty(H{cluster_iter,1}) %this means that we have less that 3 nodes in this cluster
           %make a small circle if a single point
           nodes =[];Tri=[];
            if length(x{cluster_iter,1})==1 %
                nodes(:,1) = x{cluster_iter,1}+mindist*cosd(0:45:366-45);
                nodes(:,2) = y{cluster_iter,1}+mindist*sind(0:45:366-45);
                Hnew{incr_cluster,2} = 2; %the smallest possible mesh has been created!
                %total_add_true = 1;
                Tri = delaunayTriangulation([x{cluster_iter,1};nodes(:,1)],[y{cluster_iter,1};nodes(:,2)]);
                Hnew{incr_cluster,1} = Tri.ConnectivityList;
                 xnew{incr_cluster,1} = [x{cluster_iter,1};nodes(:,1)];
                 ynew{incr_cluster,1} =  [y{cluster_iter,1};nodes(:,2)];
            else % make a small circle using the line segment
                rn = [x{cluster_iter,1}(1)-x{cluster_iter,1}(2) y{cluster_iter,1}(1)-y{cluster_iter,1}(2)];
                rnp= [-rn(2), rn(1)];
                rnp = rnp/norm(rnp);
                dd=pdist2([x{cluster_iter,1}(1) y{cluster_iter,1}(1)],[x{cluster_iter,1}(2) y{cluster_iter,1}(2)]);
                Radi=0.5*dd;
                nodes(:,1) = [mean(x{cluster_iter,1});mean(x{cluster_iter,1})+Radi*cosd([0:45:361-45])'];
                nodes(:,2) = [mean(y{cluster_iter,1});mean(y{cluster_iter,1})+Radi*sind([0:45:361-45])'];
                Tri = delaunayTriangulation(nodes(:,1),nodes(:,2));
               
                Hnew{incr_cluster,1} = Tri.ConnectivityList;
                if Radi<mindist
                    Hnew{incr_cluster,2} = 0;
                else
                    Hnew{incr_cluster,2} = 1;
                    total_add_true = 1;
                end
                xnew{incr_cluster,1} = nodes(:,1);
                ynew{incr_cluster,1} = nodes(:,2);
            end
            incr_cluster = incr_cluster+1;
            
        else
            Hnew{incr_cluster,1} = H{cluster_iter,1};
            Hnew{incr_cluster,2} = 0; %this means that it has reached the minimize size that this small region can reach!
            xnew{incr_cluster,1} = x{cluster_iter,1};
            ynew{incr_cluster,1} = y{cluster_iter,1};
            incr_cluster = incr_cluster+1;
                        
        end 
    end
        

end
ynew( ~cellfun('isempty',xnew)==0) = [];
Hnew(~cellfun('isempty',xnew)==0,:)=[];
xnew((~cellfun('isempty',xnew)==0)) = [];  


                    


%This function takes the middle point of the edges of the triangles
function  [extra_nodes,add_node] = fun_add_nodes(x,y,H,MinDist)
add_node = 0;
extra_nodes = [];
 if size(x,1)>=3
        Hsub = H;
        
   %     extra_temp(:,1) = [x1(Hsub(:,1))+0.5*(x1(Hsub(:,2))- x1(Hsub(:,1)));x1(Hsub(:,1))+0.5*(x1(Hsub(:,3))- x1(Hsub(:,1)));x1(Hsub(:,3))+0.5*(x1(Hsub(:,2))- x1(Hsub(:,3)))];
    %    extra_temp(:,2) = [y1(Hsub(:,1))+0.5*(y1(Hsub(:,2))- y1(Hsub(:,1)));y1(Hsub(:,1))+0.5*(y1(Hsub(:,3))- y1(Hsub(:,1)));y1(Hsub(:,3))+0.5*(y1(Hsub(:,2))- y1(Hsub(:,3)))];
       
        if size(Hsub,1) == 1
               extra_nodes_temp = [mean(x(Hsub)) mean(y(Hsub))];

               limD = min(pdist2(extra_nodes_temp,[x y]));
               if limD>=MinDist
                   add_node = 1;
                  extra_nodes = extra_nodes_temp;

               end
           
        else
              extra_nodes_temp = [mean(x(Hsub),2) mean(y(Hsub),2)];

              D= pdist2(extra_nodes_temp,[x y]);
              closeD =  min(D,[],2);
           %  plot( extra_nodes_temp(:,1), extra_nodes_temp(:,2),'.r')
             
             
             % limD = median(closeD);

              %Keep only the points which are above the limit limD
              %if the limD<MinDist then do not add any point

              if max(closeD)>=MinDist
                   add_node = 1;
                  extra_nodes = extra_nodes_temp(closeD>=MinDist,:);
                  
                
            

              end
        end
    
    elseif length(x)==2
        
        extra_nodes(:,1) = x(1)+0.5*(x(2)-x(1));
        extra_nodes(:,2) = y(1)+0.5*(y(2)-y(1));
    else
        extra_nodes = [];
      
        
 end 
 

 %this function divides a given mesh to independant smalles meshes based on
%the clustering of the non-zero nodes
function [Hdiv,xdiv,ydiv] = fun_div_mesh(x,y,Ind,H,cluster_number,clusters)


ind_nodes_in_a_cluster = find(abs(Ind(:,2)- cluster_number)<1e-8);      


Hsub = H(clusters(clusters(:,2)==cluster_number,1),:);
Hdiv = [];
if isempty(Hsub)
    
    xdiv = x(Ind(ind_nodes_in_a_cluster,1),1); 
    ydiv = y(Ind(ind_nodes_in_a_cluster,1),1); 
    
    
else

    % update the indices of the elements according to the created clusters
   [data,Hdiv] = RearrangeTri(1000*Hsub ,x,y);
     
    xdiv = data(:,1);
    ydiv = data(:,2);

end


function [points,Hn] = RearrangeTri(Ho,x,y)

Hn = zeros(size(Ho));
IndHo = unique(sort(Ho(:))); 
points = zeros(length(IndHo),2);
    for k = 1:length(IndHo)

       ind_temp = find(Ho(:) == IndHo(k));
       points(k,1)= x(IndHo(k)*1e-3,1);
       points(k,2)= y(IndHo(k)*1e-3,1);
       Hn(ind_temp) = k;
    end
Hn = reshape(Hn,size(Ho,1),size(Ho,2));



