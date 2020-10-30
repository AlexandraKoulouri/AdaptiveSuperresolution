function [c_es,xi_es,Ind] = fun_EstPeakLocation_noBrIt(x,c_hat,sigma,lambda,aa)
%return the amplitude c_es and the location xi_es of the estimated peak
%return also the indices of the grid points next to the peaks
%h = abs(x(2)-x(1));
N = length(x);


c_es  = [];
xi_es = [];
l_es  = [];
Ind   = [];
iter  = 1;
IndPeak = 1;
bool = 0;% bool is zero when there is only one non-zero coefficient at a grid point

while ~isempty(IndPeak)
    
    IndPeak = find(c_hat==max(c_hat) & c_hat>1e-3);
    
    if isempty(IndPeak)  break;  end
    
    for i = 1:length(IndPeak)
    
        if IndPeak == 1
            a = c_hat(1);
            
            %Check the neiboring point
            b = c_hat(2);
            c_hat(1) = 0; c_hat(2) = 0;
            if b>1e-4
                bool =1 ; %keep this point
                Indtemp =[1 2];
                xa = x(1);
                h = x(2) - x(1);
            else
                bool = 0;
                xi_es(iter) = x(IndPeak);
                c_es(iter) = a;
                Indtemp = 1;
                iter = iter+1;
            end
           
            
        elseif IndPeak == N
            a = c_hat(N-1); b = c_hat(N);
            c_hat(N-1) = 0; c_hat(N) = 0;
            
            if a>1e-4
                bool = 1; %keep the point 
                xa = x(N-1);
                Indtemp =[N-1 N];
                h = x(N) - x(N-1);
            else
                bool = 0;
                xi_es(iter) = x(N);
                c_es(iter) = b;
                Indtemp = N;
                iter = iter+1;
                
            end
            
            
        else 
            
             if c_hat(IndPeak-1)>c_hat(IndPeak+1)
                 a = c_hat(IndPeak-1);
                 b = c_hat(IndPeak);
                 c_hat(IndPeak) = 0;
                 c_hat(IndPeak-1) = 0;
                
                 if a>1e-4
                     bool = 1;
                     xa = x( IndPeak-1);
                     Indtemp =[ IndPeak-1 IndPeak];
                     h = x(IndPeak) - x(IndPeak-1);
                 else
                     bool = 0;
                     xi_es(iter) = x(IndPeak);
                     c_es(iter) = b;
                     Indtemp =[IndPeak];
                     iter = iter+1;
                 end
             else
                b = c_hat(IndPeak+1);
                a = c_hat(IndPeak);
                c_hat(IndPeak) = 0;
                c_hat(IndPeak+1) = 0;
                if b>1e-4              
                   bool = 1;
                   Indtemp =[ IndPeak IndPeak+1];
                   xa = x( IndPeak);
                   h = x(IndPeak+1) - x(IndPeak);
                else
                    bool = 0;
                    xi_es(iter) = x(IndPeak);
                    c_es(iter) = a;
                    iter = iter+1;
                end
            end
        end
        
        if bool == 1
            [xi_es(iter), c_es(iter)] = EstimateLocation (xa, a,b,sigma,h,lambda,aa); %#ok<AGROW>
             iter = iter + 1;
             bool = 0;
        end
        Ind = [Ind Indtemp];
     end
        
   
end



function [xi_es, c_es] = EstimateLocation (xa, a,b,sigma,h,lambda,aa) %#ok<DEFNU>
% sigma = sqrt(2)*sigma;
% % G_double = @(u) aa*2*sigma*sqrt(pi)*normpdf(u,0,sigma_double);
%         
% 
 G_double = @(u) aa*exp(-0.5 * ( u./sigma).^2)./ (sqrt(2*pi) .* sigma);
c_es = a+b+lambda/(2*G_double(0));
xi_es = 0.5*(2*xa+h) + (b-a)*h/(2*c_es); 

%xi_es = 0.5*(h*(a-b)/(a+b)+2*xa+h );
%c_es = a+b; 

%G_double = @(u) aa*exp(-0.5 * ( u./sigma).^2)./ (sqrt(2*pi) .* sigma);

%la_es = 0.5- ((a-b)/(a+b)*(G_double(0)-G_double(h)))/(h^2*G_double(0)/sigma^2);


%xi_es = xa+ la_es*h;
%c_es = a+b;

% xi_es = 0.5*(h*(a-b)/c_es+2*xa+h );
%c_es = a+b; 
 
