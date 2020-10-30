function A = STORM_2D_Gen_Meas_Mat_v1(xxx_r,yyy_r,xx_r,yy_r,sigma)


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