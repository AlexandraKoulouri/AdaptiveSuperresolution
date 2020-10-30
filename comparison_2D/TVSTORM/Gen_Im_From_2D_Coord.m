function im = Gen_Im_From_2D_Coord(coord_list, I_list, height, width, sigma)
% PURPOSE:
% Generate image based on input emitter coordinate and intensity.
%---------------------------------------------------
% USAGE:
% im = Gen_Im_From_Coord(coord_list, I_list, height, width, sigma)
%---------------------------------------------------
% INPUTS:
% coord_list:   list of input coordinate
% I_list:       list of input intensity
% height:       image height
% width:        image width
% sigma:        standard deviation of system PSF (pixel)
%---------------------------------------------------
% OUTPUTS:
% im:           generated image
%---------------------------------------------------

[xx, yy] = meshgrid(1 : width, 1 : height);
xx_r = xx(:); yy_r = yy(:);

im = zeros(height * width, 1);

% iterate through all input coordinate
for i = 1 : size(coord_list, 1)
    % generate measurement on a given point
    im = im + I_list(i) / 4 *...
        (erf((xx_r + 0.5 - coord_list(i, 1)) / sqrt(2) / sigma) -...
        erf((xx_r - 0.5 - coord_list(i, 1)) / sqrt(2) / sigma)) .*...
        (erf((yy_r + 0.5 - coord_list(i, 2)) / sqrt(2) / sigma) -...
        erf((yy_r - 0.5 - coord_list(i, 2)) / sqrt(2) / sigma));
end