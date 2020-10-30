function [x,y,Ind_nz_cell,c] = removeEmptyCells(x,y,Ind_nz_cell,c);
%Ind_empty = (~cellfun('isempty',c_hat_cell)==0);
y( ~cellfun('isempty',c)==0) = [];
x( find(~cellfun('isempty',c)==0)) = [];
Ind_nz_cell(find(~cellfun('isempty',c)==0))=[];
c(find(~cellfun('isempty',c)==0))=[];

y (find(~cellfun('isempty',x)==0)) = [];
c(find(~cellfun('isempty',x)==0)) = [];
Ind_nz_cell(find(~cellfun('isempty',x)==0))=[];
x (~cellfun('isempty',x)==0) = [];   