function A = GenerateBasisFunMat2D(Gpoints,x,sigma)



%Basis functions (constant sigma)

Num1 = size(Gpoints,1);
Num2 = size(x,1); %canditate grid points

sigma_x = sigma;
sigma_y = sigma;


%G = @(u1,u2)1/(2*pi*(sigma_x*sigma_y))*exp(-0.5*(u1.^2/sigma_x^2+u2.^2/sigma_y^2));
%G = 


u1 = repmat(x(:,1)',Num1,1) - repmat(Gpoints(:,1),1,Num2);
u2 = repmat(x(:,2)',Num1,1) - repmat(Gpoints(:,2),1,Num2);
A = GaussianKernel2D(sigma_x,sigma_y,u1,u2);


%Basis functions (constant sigma)

% Num1 = size(Gpoints,1);
% Num2 = size(x,1); %canditate grid points
% 
% sigma_x = sigma;
% sigma_y = sigma;
% theta = 0;
% 
% a = cos(theta)^2/2/sigma_x^2 + sin(theta)^2/2/sigma_y^2;
% b = -sin(2*theta)/4/sigma_x^2 + sin(2*theta)/4/sigma_y^2 ;
% c = sin(theta)^2/2/sigma_x^2 + cos(theta)^2/2/sigma_y^2;
% 
% G = @(u1,u2)1/(2*pi*(sigma_x*sigma_y))* exp( - (a*(u1).^2 - 2*b*(u1).*(u2) + c*(u2).^2)) ;
% 
% 
% 
% u1 = repmat(x(:,1)',Num1,1) - repmat(Gpoints(:,1),1,Num2);
% %u2 = repmat(y0,1,length(y0)) - repmat(y0',length(y0),1);
% u2 = repmat(x(:,2)',Num1,1) - repmat(Gpoints(:,2),1,Num2);
% A = G(u1, u2);