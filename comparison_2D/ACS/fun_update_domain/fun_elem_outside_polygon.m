function [t_rem] =fun_elem_outside_polygon(g,mid_p, t_ind)


%remove elements
ind_out = find(inpolygon(mid_p(:,1),mid_p(:,2),g(:,1),g(:,2))==0);

t_rem = unique(t_ind(ind_out));
end