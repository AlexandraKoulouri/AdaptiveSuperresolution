function [clustCent, xm] = STORM_2D_CS(A, im, height, width, div,...
    thresh, ns, di_wid)
% PURPOSE:
% STORM 2D reconstruction using compressive sensing.
%---------------------------------------------------
% USAGE:
% [clustCent, xm] = STORM_2D_CS(A, im, height, width, div,...
%   thresh, ns, di_wid)
%---------------------------------------------------
% INPUTS:
% A:                measurement matrix
% im:               acquired image
% height:           image height
% width:            image width
% div:              upsampling factor in x y-dimension
% thresh:           threshold for reconstructed image
% ns:               parameter that controls the balance between
%                   data fidelity and sparsity
%                   typically 1.5 is used for Poisson noise
% di_wid:           dilation width in x y-dimension
%---------------------------------------------------
% OUTPUTS:
% clustCent(:, 1):  x coordinates (pixel)
% clustCent(:, 2):  y coordinates (pixel)
% xm:               reconstructed upsampled image
%---------------------------------------------------

%% Default values
% if nargin <= 9
%     ns = 1.5;
% end
% if nargin <= 10
%     di_wid = 5;
% end

%% Pre-processing for compressed sensing
im = double(im(:));
eps = ns * sqrt(sum(im));

%% Compressed sensing on the image
xm = STORM_Homotopy(A, im, 'tol', eps,'isnonnegative', true, 'maxiteration', 500);

%% Process the result
xm = xm(1 : height * div * width * div);
x2 = reshape(xm, height * div, width * div);
% thresholding
BW = x2 > thresh;
se = ones([di_wid, di_wid]);
% dilation
BW = imdilate(BW, se);
% debiasing
idx = find(BW(:) ~= 0);
A2 = A(:, idx);
xmm = xm;
xmm(idx) = lsqnonneg(A2, im);
% weighted centroid
x2 = reshape(xmm, height * div, width * div);
STAT = regionprops(BW == 1, x2, 'WeightedCentroid', 'PixelValues');

clustCent = cat(1, STAT.WeightedCentroid);
for statid = 1 : size(STAT, 1)
    clustCent(statid, 3) = sum(STAT(statid).PixelValues);
end

if ~isempty(clustCent)
    clustCent(isnan(clustCent(:, 1)), :) = []; % remove nan
    clustCent(:, 1) = (clustCent(:, 1) + (div / 2 - 0.5)) / div;
    clustCent(:, 2) = (clustCent(:, 2) + (div / 2 - 0.5)) / div;
end