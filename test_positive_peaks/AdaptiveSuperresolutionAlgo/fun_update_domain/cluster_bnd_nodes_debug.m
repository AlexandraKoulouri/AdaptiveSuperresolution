function [Ind_Edge,e] = cluster_bnd_nodes_debug(e,GroupNum,x,y)
% e : are the indices of the nodes along the edges
%Groupd Num: the group in which the nodes belongs 1,2,.....
[foo,ix,jx] = unique(e(:));
vec = histc(jx,1:max(jx)); 
xx  = find(vec == 4); 

%e= sort(e,2);

% if there are two regions  which are connected with a single node or with
% an element(triangle)/small region
if ~isempty(xx) && length(xx) >1 
    t = foo(xx(1));
    xx(1) = [];
    indt = find(e(:,1)==t(1)); %check both columns to find the node
    if isempty(indt)
       indt = find(e(:,2)==t(1));
    end
       %first node
       
       Ind_Edge(1,1) = e(indt(1),1);
       Ind_Edge(2,1) = e(indt(1),2);
       Ind_Edge(1,2) = GroupNum; 
       Ind_Edge(2,2) = GroupNum; 
            
     % Ext_Edge = find( xx == Ind_Edge(2,1));
      
       if ~isempty(find(foo(xx) == Ind_Edge(2,1)))
           CurrentEdgeP = Ind_Edge(1,1); 
           Term_Edge    = Ind_Edge(2,1);
           Ind_intersection =3;
       else          
          CurrentEdgeP  = Ind_Edge(2,1);
          Term_Edge     = Ind_Edge(1,1);
       end
      e(indt(1),:) = [];
elseif( ~isempty(xx) && length(xx) ==1 )%if there is only a single connection point between two regions
    t = foo(xx(1));
    xx(1) = [];
    indt = find(e(:,1)==t(1)); %check both columns to find the node
    if isempty(indt)
          indt = find(e(:,2)==t(1));
    end


%     CurrentEdgeP = e(indt(1),2);
    Ind_Edge(1,1) = e(indt(1),1);
    Ind_Edge(2,1) = e(indt(1),2);
    Ind_Edge(1,2) = GroupNum; 
    Ind_Edge(2,2) = GroupNum;
    CurrentEdgeP  = Ind_Edge(2,1);
    Term_Edge     = Ind_Edge(1,1);
    
    e(indt(1),:) = [];
    % elseif  ~isempty(xx) && length(xx) == 2 % if there is an edge that
    % separates the two regions

else
    %CurrentEdgeP = e(1,2);
    Ind_Edge(1,1) = e(1,1);
    Ind_Edge(2,1) = e(1,2);
    Ind_Edge(1,2) = GroupNum; 
    Ind_Edge(2,2) = GroupNum;
    CurrentEdgeP  = Ind_Edge(2,1);
    Term_Edge     = Ind_Edge(1,1);
    
    e(1,:) = [];

end
iter = 3; 
%hold on
%plot(x(Ind_Edge(:,1)),y(Ind_Edge(:,1)),'r')
while ~isempty(e)  
 
    IndnextEdge =  find(e(:,1) == CurrentEdgeP);
    if ~isempty(IndnextEdge) 
           if(e(IndnextEdge(1),2) ~= Term_Edge)
               Ind_Edge(iter,1) = e(IndnextEdge(1),2);
               Ind_Edge(iter,2) = GroupNum;
               CurrentEdgeP = e(IndnextEdge(1),2);
               iter = iter+1; 
               if ~isempty(xx)
                    if ~isempty(find(foo(xx)==CurrentEdgeP))
                        Ind_Term_Edge    = find(foo(xx) == CurrentEdgeP);
                        Term_Edge  = foo(xx(Ind_Term_Edge));
                        CurrentEdgeP = Ind_Edge(1,1);
                        Ind_intersection = iter;
                    elseif (CurrentEdgeP==Ind_Edge(1,1))
                        %Ind_brGroups = find(Ind_Edge(:,1) == CurrentEdgeP); % divide into two groups
                        Ind_Edge(Ind_intersection:end,2) = GroupNum +1;  %change group number
                        Ind_intersection=[];
                        CurrentEdgeP = Ind_Edge(1,1); %go back to the first group
                        %GroupNum = GroupNum -1;
                        
                    end              
               end
                
            elseif (e(IndnextEdge(1),2) == Term_Edge || Term_Edge~=Ind_Edge(1,1) )
               % Ind_Edge(iter,1) = e(IndnextEdge(1),2);
               % Ind_Edge(iter,2) = GroupNum;
                iter = iter+1;
            end
           e(IndnextEdge(1),:) = [];
           
   
    else
        IndnextEdge =  find(e(:,2) == CurrentEdgeP); 
         if ~isempty(IndnextEdge) 
            %swap the nodes of the edge
            temp = e(IndnextEdge(1),2);
            e(IndnextEdge(1),2) = e(IndnextEdge(1),1);
            e(IndnextEdge(1),1) = temp;
         else
            %create a new group
            GroupNum = max(Ind_Edge(:,2))+1;%GroupNum +1;
            [Ind_Group_t,e] = cluster_bnd_nodes_debug(e,GroupNum,x,y);
            Ind_Edge =[Ind_Edge;Ind_Group_t]; %#ok<AGROW>
         end  
    end
    %   plot(x(Ind_Edge(1:26,1)),y(Ind_Edge(1:26,1)),'g')
    

end


