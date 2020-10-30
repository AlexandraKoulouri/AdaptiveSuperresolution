function A = DirichletKernel(fc,x,z)

z_rep = repmat(z,1,length(x));
x_rep = repmat(x,length(z),1);

A = (sin((2*fc+1)*pi*(z_rep-x_rep))./sin(pi*(z_rep-x_rep)))/(2*fc+1);

A(isnan(A))=1;