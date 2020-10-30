function A = STORM_2D_Gen_Meas_Mat(width, height, div, sigma)
% PURPOSE:
% This function is to generate the measurement matrix for 2D single
% camera STORM image reconstruction.
%---------------------------------------------------
% USAGE:
% A = STORM_2D_Gen_Meas_Mat(width, height, div, sigma)
%---------------------------------------------------
% INPUTS:
% width:        image width
% height:       image height
% div:          upsampling factor in x y-dimension
% sigma:        standard deviation of system PSF (pixel)
%---------------------------------------------------
% OUTPUTS:
% A:            generated measurement matrix
%---------------------------------------------------

% meshgrid of original image size
[xx, yy] = meshgrid(1 : width, 1 : height);
% meshgrid of upsampled image size
[xxx, yyy] = meshgrid(0.5 + 1 / div / 2 : 1 / div : width + 0.5 - 1 / div / 2,...
    0.5 + 1 / div / 2 : 1 / div : height + 0.5 - 1 / div / 2);
xx_r = xx(:); yy_r = yy(:);
xxx_r = xxx(:); yyy_r = yyy(:);

A = zeros(width * height, size(xxx_r, 1));

% iterate through all upsampled grid points
for i = 1 : size(xxx_r, 1)
    cx = xxx_r(i);
    cy = yyy_r(i);
    % generate measurement on a given point
    A(:,i) = 1 / 4 *...
        (erf((xx_r + 0.5 - cx) / sqrt(2) / sigma) -...
        erf((xx_r - 0.5 - cx) / sqrt(2) / sigma)) .*...
        (erf((yy_r + 0.5 - cy) / sqrt(2) / sigma) -...
        erf((yy_r - 0.5 - cy) / sqrt(2) / sigma));
end