function [y,xi,h0] = SimulateObservations(x,sigma,mu_i,c)

%Generate gaussian convoluted data
%G = @(u) sqrt(2*pi)*sigma*normpdf(u,0,sigma);
G = @(u)  exp(-0.5 * ((u)./sigma).^2) ./ (sqrt(2*pi) .* sigma);


%amplitudes

%Location of the peak
xi = mu_i;
k = length(mu_i);
y = zeros(length(x),1);
h = min(diff(x));
h0 = h(1);

for i = 1:k
        %Observation
        y = y+ c(i)*G(x-mu_i(i))';
end
        



   







