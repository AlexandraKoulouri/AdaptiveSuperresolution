function [t_rem] =fun_elem_outside_polygon(p,mid_p, t_ind,tri)


%remove elements
ind_out = find(inpolygon(mid_p(:,1),mid_p(:,2),p(:,1),p(:,2))==0);

t_rem = unique(t_ind(ind_out));

%check if the element (about to remove) have a node which does not belong
%to any other element (then keep this element)



if ~isempty(t_rem)
    ind_t_left = (setdiff(1:size(tri,1),t_rem));
    t_left = tri(ind_t_left,:);
    keep_t = [];
    for i=1:size(t_rem,1)

        nodes_temp_ind = tri(t_rem(i),:);

       test1=find(t_left(:,1)==nodes_temp_ind (1) |  t_left(:,2)==nodes_temp_ind (1) |  t_left(:,3)==nodes_temp_ind (1));
       test2=find(t_left(:,1)==nodes_temp_ind (2) |  t_left(:,2)==nodes_temp_ind (2) |  t_left(:,3)==nodes_temp_ind (2));
       test3=find(t_left(:,1)==nodes_temp_ind (3) |  t_left(:,2)==nodes_temp_ind (3) |  t_left(:,3)==nodes_temp_ind (3));
       if isempty(test1) ||  isempty(test2) || isempty(test3)
           keep_t = [keep_t;i];
       end


    end
t_rem(keep_t)=[];
end

