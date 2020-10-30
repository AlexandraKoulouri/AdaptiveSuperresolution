function ob= SimulateObservationsGaussian(x,y,sig,amp,mu_xi,mu_yi,alpha);

%Generate gaussian convoluted data give us
% y= (a1*G(x-mu1)+a2*G(x-mu1)+...)+(a1*G(x-mu2)+a2*G(x-mu2)+...)

ob = zeros(length(x),1);
%Gaussian kernel
G = @(u1,u2,sigma)1/(2*pi*(sigma^2))*exp(-0.5*(u1.^2/sigma^2+u2.^2/sigma^2));

for i = 1:length(sig)
    for j = 1:length(mu_xi);
      ob = ob+alpha(i)*amp(j)*G(x-mu_xi(j),y-mu_yi(j),sig(i));
    end
end
ob=ob';

   
% function img_kernel = Gauss_kernel(x_inx,y_inx,x_pos,y_pos,Gsigma)
%     img_kernel = exp(-((x_inx-x_pos).^2+(y_inx-y_pos).^2)/(2*Gsigma^2))/(2*pi*Gsigma^2) ;
% end