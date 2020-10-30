function [xi_es,c_approx] = fun_loc_est(c,x,y,Sigma,alpha,lambda)
%2D cases

%c    : coefficients
%(x,y): coordinates of the non-zero coefficients
%returns the locations of the peak

%Check collinear points first!If this is the case then solve the problem in
%one dimension along the line (this can be evolved by solving the problem
%in one dimension but here we will do it very simply
iter1 = 1;
if(length(x)>=3)
    tf = rank(bsxfun(@minus, [x y], [x(1,:) y(1,:)]), eps) < 2;
    ctemp = c;
    if tf == 1 %then the points are collinear
        ind = zeros(length(x)-2,1);
        for iter = length(x):-1:3
            ind(iter1) = find(ctemp==min(ctemp));
            ctemp(ind(iter1)) = max(ctemp);
            iter1 = iter1+1;
        end

        %remove the point locations with the small amplitude values

    x(ind)=[];
    y(ind)=[];
    end
end



  NumOfTerms = length(alpha);
    Grad2_H_0 = zeros(NumOfTerms);
    H_0=[];
for i = 1:length(alpha)
    
      [SecD, H_0_t] = Gaussian_2Der(sqrt(2)*Sigma(i),0,0);
        
     Grad2_H_0 = Grad2_H_0 + alpha(i)^2*SecD;
     H_0 = H_0+H_0_t;
  
    
end
if NumOfTerms>1
    Ind_a_comb = combnk(1:NumOfTerms,2);
    for i = 1:size(Ind_a_comb,1)
        Grad2_H_0 = Grad2_H_0+2*alpha(Ind_a_comb(i,1))*alpha(Ind_a_comb(i,2))*Gaussian_2Der(sqrt(Sigma(Ind_a_comb(i,1))+Sigma(Ind_a_comb(i,2))),0,0);
    end
    
end

%approximate amplitude
c_approx = sum(c);%+lambda/(2*H_0);
xi_es = zeros(2,1);
P = [x y];

%include the case of co-linearity (I haven't consider it yet)
if length(x)>=3

    ind = combnk(1:length(x),2); 

    
    b = zeros(length(ind),1);
    for i =1:length(ind);

        b(i) = 0.5 * ( P(ind(i,1),:)*Grad2_H_0*P(ind(i,1),:)' - P(ind(i,2),:)*Grad2_H_0*P(ind(i,2),:)');

        for k = 1:length(x)

             b(i) = b(i)-0.5/c_approx*c(k)*( (P(ind(i,1),:)-P(k,:))*Grad2_H_0*(P(ind(i,1),:)-P(k,:))' -...
                 (P(ind(i,2),:)-P(k,:))*Grad2_H_0*(P(ind(i,2),:)-P(k,:))');
            %  b(i) = b(i) + c(k)*P(k,:)*Grad2_H_0*(P(ind(i,1),:)-P(ind(i,2),:))';
        end

    end 
  
    A = zeros(2,length(ind)); %#ok<PREALL>
    A = (P(ind(:,1),:)-P(ind(:,2),:))*Grad2_H_0;
    xi_es = A\b ;
elseif length(x) == 2
    %estimate the location along a line
    %parametrize the the xi with a line (xi+t(xj-xi))%Maybe there is an
    %error
      b =0;
      b = P(1,:)*Grad2_H_0*P(1,:)'-P(2,:)*Grad2_H_0*P(2,:)';
      b = 0.5*b - 0.5/c_approx*(-c(1)*(P(2,:)-P(1,:)) *Grad2_H_0*(P(2,:)-P(1,:))'+c(2)*(P(1,:)-P(2,:))*Grad2_H_0*(P(1,:)-P(2,:))');

     Aa = (P(1,:)-P(2,:))*Grad2_H_0;
     t0=  Aa*(P(1,:)-P(2,:))';
     t1= Aa*P(2,:)';
     t = (b-t1)/t0;        
              
     xi_es = P(2,:)+t*(P(1,:)-P(2,:));
    
else
   %if only one non-zero point then assign the closest grid point location
     xi_es = [x y]';
    
end

function [G_2,G_0] = Gaussian_2Der(sigma,x,y)
%non-correlation between sigma_x and sigma_y

sigma_x = sigma;
sigma_y = sigma;

% %Gaussian kernel 
G = @(u1,u2)1/(2*pi*(sigma_x*sigma_y))*exp(-0.5*(u1.^2/sigma_x^2+u2.^2/sigma_y^2));
G_0 = G(x,y);

%first derivative
% G_p_x = @(u1,u2)-1/(sigma_x^2*2*pi*(sigma_x*sigma_y))*exp(-0.5*(u1.^2/sigma_x^2+u2.^2/sigma_y^2)).*u1;
% G_p_y = @(u1,u2)-1/(2*pi*(sigma_y^2*sigma_x*sigma_y))*exp(-0.5*(u1.^2/sigma_x^2+u2.^2/sigma_y^2)).*u2;
% G_p = [G_p_x(x,y), G_p_y(x,y)];

%Second derivative
G_2_x = @(u1,u2)-1/(sigma_x^2*2*pi*(sigma_x*sigma_y))*exp(-(u1.^2/sigma_x^2+u2.^2/sigma_y^2))...
    +2/(sigma_x^4*2*pi*(sigma_x*sigma_y))*exp(-(u1.^2/sigma_x^2+u2.^2/sigma_y^2))*u1;
G_2_y =  @(u1,u2)-1/(sigma_y^2*2*pi*(sigma_x*sigma_y))*exp(-(u1.^2/sigma_x^2+u2.^2/sigma_y^2))...
    +2/(sigma_y^4*2*pi*(sigma_x*sigma_y))*exp(-(u1.^2/sigma_x^2+u2.^2/sigma_y^2))*u2;

G_2_xy= @(u1,u2)-1/(sigma_x^2*sigma_y^2*2*pi*(sigma_x*sigma_y))*exp(-(u1.^2/sigma_x^2+u2.^2/sigma_y^2)).*u1.*u2;

G_2 = [G_2_x(x,y) G_2_xy(x,y); G_2_xy(x,y)' G_2_y(x,y)];
    
%  function [h,s,A] = Triangle(x,y)
% s = zeros(3,1);
% s(1) = sqrt((x(1)-x(2)).^2+(y(1)-y(2)).^2);
% s(2) = sqrt((x(2)-x(3)).^2+(y(2)-y(3)).^2);
% s(3) = sqrt((x(1)-x(3)).^2+(y(1)-y(3)).^2);
% p = 1/2*sum(s); %triangle surface
% A =sqrt(p*(p-s(1))*(p-s(2))*(p-s(3)));
% h = 2*A./s;%triangle hights
    