%% ****************************************************************************************************
%% each localized particle is matched to the closest true particle within a radius
%% ****************************************************************************************************
function [Recall1,Precision1,F_score1,Jaccard1] = EvaluationMetrics_scores(Results,true_poses,true_amplitude,h,radius)

% the domain is [0,1]^2
% h is defined by 1/size of image If not insert then just but h = 1;

TP = 0;              %number of correctly localized
FP = 0;
FN = 0;
est_pos = Results(:,2:3);                  %Estimated positions
true_pos = true_poses(:,:);             %Actual positions
Na = size(est_pos,1);                   %number of reconstructed peaks  
Nb = size(true_amplitude,1);            %number of original peaks
%find minimum distance between estimated positions

if radius == 0
    dist = zeros(size(est_pos,1));
    for i = 1:size(est_pos,1)-1;
       for j = i+1:size(est_pos,1)
          dist(i,j)=sqrt((est_pos(i,1)-est_pos(j,1)).^2 + (est_pos(i,2)-est_pos(j,2)).^2);
          dist(j,i) = dist(i,j);
       end

    end
    mdist = min(dist(:));
    radius = 0.01*mdist; %tolerance
end


indInsideCircle = []; %#ok<NASGU>
for ii = 1:size(est_pos,1)
        % find the closest true location
        distance = sqrt((true_pos(:,1)-est_pos(ii,1)).^2 + (true_pos(:,2)-est_pos(ii,2)).^2); %find the minimum distance
        indInsideCircle = find(distance<=radius);
        if isempty(indInsideCircle)
            FP = FP+1;
        else
            if length(indInsideCircle)==1
              TP = TP+1;
              true_pos(indInsideCircle,:) = []; %remove the true point
            else
               TP = TP+1;
               FN = FN+length(indInsideCircle)-1;
               true_pos(indInsideCircle,:) = [];
                
            end
            
        indInsideCircle = [];
      end
    
  end
if ~isempty(true_pos)
    FN = FN + length(true_pos,1);
end
%Metrics
Recall1 = TP/(TP+FN);
Precision1  = TP/(TP+FP);
F_score1 = 2*(Recall1*Precision1)/(TP+FN);
Jaccard1 = TP/(Na+Nb-TP);
