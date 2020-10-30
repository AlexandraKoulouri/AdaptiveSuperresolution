function A = Gen_Meas_From_2D_Coord(coord_list, height, width, sigma)
% PURPOSE:
% Generate measurement matrix based on input emitter coordinate.
%---------------------------------------------------
% USAGE:
% A = Gen_Meas_From_2D_Coord(coord_list, height, width, sigma)
%---------------------------------------------------
% INPUTS:
% coord_list:   list of input coordinate
% height:       image height
% width:        image width
% sigma:        standard deviation of system PSF (pixel)
%---------------------------------------------------
% OUTPUTS:
% A:            generated measurement matrix
%---------------------------------------------------

[xx, yy] = meshgrid(1 : width, 1 : height);
xx_r = xx(:); yy_r = yy(:);

A = zeros(height * width, size(coord_list, 1));

% iterate through all input coordinate
for i = 1 : size(coord_list, 1)
    % generate measurement on a given point
    A(:, i) = 1 / 4 *...
        (erf((xx_r + 0.5 - coord_list(i, 1)) / sqrt(2) / sigma) -...
        erf((xx_r - 0.5 - coord_list(i, 1)) / sqrt(2) / sigma)) .*...
        (erf((yy_r + 0.5 - coord_list(i, 2)) / sqrt(2) / sigma) -...
        erf((yy_r - 0.5 - coord_list(i, 2)) / sqrt(2) / sigma));
end