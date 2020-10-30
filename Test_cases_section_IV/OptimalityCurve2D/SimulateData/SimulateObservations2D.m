function ob= SimulateObservations2D(x,y,sigma,amp,mu_xi,mu_yi);

%Generate gaussian convoluted data
theta = 0;

sigma_x = sigma;
sigma_y = sigma;

a = cos(theta)^2/2/sigma_x^2 + sin(theta)^2/2/sigma_y^2;
b = -sin(2*theta)/4/sigma_x^2 + sin(2*theta)/4/sigma_y^2 ;
c = sin(theta)^2/2/sigma_x^2 + cos(theta)^2/2/sigma_y^2;

%G = @(u1,u2) 1/(2*pi*(sigma_x*sigma_y))* exp( - (a*(u1).^2 - 2*b*(u1).*(u2) + c*(u2).^2)) ;
G = @(u1,u2) exp( - (a*(u1).^2 - 2*b*(u1).*(u2) + c*(u2).^2)) ;

% % %Observation
ob = zeros(length(x),1);
for i = 1:length(mu_xi);
ob = ob+amp(i)*G(x-mu_xi(i),y-mu_yi(i));
end

ob=ob';

   
