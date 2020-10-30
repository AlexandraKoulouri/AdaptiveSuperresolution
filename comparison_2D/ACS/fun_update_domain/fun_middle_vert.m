function [mid_p,t_ind] = fun_middle_vert(p,t)
num_t = size(t,1);

ver   =[t(:,[1,2]);
       t(:,[1,3]);
       t(:,[2,3])];
 
   
   
 t_ind = repmat(1:num_t,1,size(t,2));  
 mid_p = (p(ver(:,1),:)-p(ver(:,2),:))/2+p(ver(:,2),:);
end
 

 
 
   
   
   