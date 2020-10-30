function [clustCent, I_list] = TVSTORM_2D(A, im, height, width, sigma,...
    div, ns, thresh)
% PURPOSE:
% STORM 2D reconstruction using TVSTORM.
%---------------------------------------------------
% USAGE:
% [clustCent, I_list] = TVSTORM(A, im, height, width, sigma, div, ns, thresh)
%---------------------------------------------------
% INPUTS:
% A:                measurement matrix
% im:               acquired image
% height:           image height
% width:            image width
% sigma:            standard deviation of system PSF (pixel)
% div:              upsampling factor in x y-dimension
% ns:               parameter that controls the balance between
%                   data fidelity and sparsity
%                   typically 1.5 is used for Poisson noise
% thresh:           threshold for backward step
%---------------------------------------------------
% OUTPUTS:
% clustCent(:, 1):  x coordinates (pixel)
% clustCent(:, 2):  y coordinates (pixel)
% I_list:           intensity list
%---------------------------------------------------

%% Pre-processing for TVSTORM
im = double(im(:));
eps = ns * sqrt(sum(im));
[xxx, yyy] = meshgrid(0.5 + 1 / div / 2 : 1 / div : width + 0.5 - 1 / div / 2,...
    0.5 + 1 / div / 2 : 1 / div : height + 0.5 - 1 / div / 2);

%% Iterate
coord_list = [];
I_list = [];
bg = min(im);
iter = 1;
max_iter = 40;
res_iter = [];

coord_list_before = [];
I_list_before = [];

while iter < max_iter
    % calculate the residual
    im_gen = Gen_Im_From_2D_Coord(coord_list, I_list, height, width, sigma);
    res = im - im_gen;
    
    res_iter(iter) = norm(res);
    % terminate if residual small enough
    if res_iter(iter) <= eps
        break;
    end
    
    if sum(I_list < thresh) ~= 0
        if sum(I_list < thresh) == 1 && I_list(end) < thresh
            coord_list = coord_list_before;
            I_list = I_list_before;
            break;
        else
            while sum(I_list < thresh) ~= 0
                list = I_list > thresh;
                coord_list = coord_list(list, :);
                I_list = I_list(list, :);
                [coord_list, I_list, bg] = Coordinate_Descent_2D(...
                    coord_list, I_list, bg, im, height, width, sigma);
            end
            break;
        end
    end
    
    coord_list_before = coord_list;
    I_list_before = I_list;
    
    % find one atom
    coeff = A' * (res - bg);
    [~, idx] = max(coeff);
    coord_list = [coord_list; xxx(idx), yyy(idx)];
    I_list = [I_list; max(coeff(idx) / norm(A(:, idx)) ^ 2, 100)];
    % coordinate descent on the support
    [coord_list, I_list, bg] = Coordinate_Descent_2D(...
        coord_list, I_list, bg, im, height, width, sigma);
    iter = iter + 1;
end
clustCent = coord_list;