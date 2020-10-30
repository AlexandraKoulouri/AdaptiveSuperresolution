function A = GenerateGaussianMat(ObPoints,ComPoints,sig,alpha)

%         sigma1, and 2     : width of Gaussian fuctions for PSF
%         sigma_ratio       : ratio for two Gaussian functions for PSF: PSF = simga_ratio*G_sigma1 + (1-simga_ratio)*G_sigma2
%         

%Basis functions (constant sigma)

Num1 = size(ObPoints,1);      % Measurements points
Num2 = size(ComPoints,1);     % nodes to estimate peaks

%Gaussian kernel
G = @(u1,u2,sigma)1/(2*pi*(sigma^2))*exp(-0.5*(u1.^2/sigma^2+u2.^2/sigma^2));

u1 = repmat(ComPoints(:,1)',Num1,1) - repmat(ObPoints(:,1),1,Num2);
u2 = repmat(ComPoints(:,2)',Num1,1) - repmat(ObPoints(:,2),1,Num2);

A = zeros(Num1,Num2);
for i= 1:length(sig)
     A = A+alpha(i)*G(u1, u2,sig(i)); 
end
%A = A./sum(A(:));