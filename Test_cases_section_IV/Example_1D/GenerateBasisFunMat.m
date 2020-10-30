function A = GenerateBasisFunMat(Gpoints,x,sigma)

%Basis functions (constant sigma)

Num1 = length(Gpoints);
Num2 = length(x);

G = @(u) exp(-0.5 * ((u)./sigma).^2)./ (sqrt(2*pi) .* sigma); %normpdf(u,0,sigma);
u = repmat(x,Num1,1) - repmat(Gpoints',1,Num2);
A = G(u);
