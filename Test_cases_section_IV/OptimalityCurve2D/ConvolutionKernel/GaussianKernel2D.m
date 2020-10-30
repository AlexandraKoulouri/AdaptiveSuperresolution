function G = GaussianKernel2D(sigma_x,sigma_y,u1,u2)
%G = 1/(2*pi*(sigma_x*sigma_y))*exp(-0.5*(u1.^2/sigma_x^2+u2.^2/sigma_y^2));
G = exp(-0.5*(u1.^2/sigma_x^2+u2.^2/sigma_y^2));