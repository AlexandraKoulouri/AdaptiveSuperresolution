%% ****************************************************************************************************
%% calculation of localization error in x,y direction.
%% each localized particle is matched to the closest true particle within a radius
%% ****************************************************************************************************
function [num_idens,num_clusters,LE,SE,DNP,Loc_x,Loc_y,Amp] = simul_eval(Results,true_poses,true_amplitude,h,radius)
% the domain is [0,1]^2
% h is defined by 1/size of image 

num_idens = zeros(1,1);              %number of correctly localized
num_clusters = zeros(1,1);           %number of estimated peaks
% errors_x = [];
% errors_y = [];

%ind = find(Results(:,1) == tt);
est_pos = Results(:,2:3);                  %Estimated positions
true_pos = true_poses(:,:);                %Actual positions
num_cluster = size(est_pos,1);             %number of peaks  or clusters

LE = zeros(size(true_poses,1),1);         %localization error
SE = zeros(size(LE));                     %Strength error
DNP = zeros(1,1);
Loc_x = zeros(size(LE));
Loc_y = zeros(size(LE));
Amp   = zeros(size(LE));
if num_cluster
    num_iden_temp = zeros(1,size(true_pos,1));
    num_cluster_temp = zeros(num_cluster,1);
    for ii = 1: num_cluster
        % find the closest true location
        distance = sqrt((true_pos(:,1)-est_pos(ii,1)).^2 + (true_pos(:,2)-est_pos(ii,2)).^2);
        [error_2D,index] = min(distance*h);
        SE(index) = SE(index)+ abs(true_amplitude(index)-Results(ii,1));
        LE(index) = LE(index)+error_2D;
        Amp(index)= mean(Amp(index)+Results(ii,1));
        Loc_x(index)  = mean(Loc_x(index)+est_pos(ii,1));
        Loc_y(index)  = mean(Loc_y(index)+est_pos(ii,2));
        if abs(error_2D) < radius
            num_iden_temp(index) = 1;
            num_cluster_temp(ii) = 1;
        end
    end
    ind = find(num_cluster_temp>0);
    num_idens = sum(num_iden_temp(:));
    num_clusters = num_cluster;

end
DNP = abs(num_clusters-size(true_pos,1)); %difference in the source reconstructions

