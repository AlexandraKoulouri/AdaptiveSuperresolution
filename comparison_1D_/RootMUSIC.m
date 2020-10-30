%use ROOT-MUSIC 
function [xi_es,c_es] = RootMUSIC(y,N,fc)

%Based on web-page: https://www.numerical-tours.com/matlab/sparsity_9b_music/
%Publication: W. Liao, A. Fannjiang, MUSIC for Single-Snapshot Spectral Estimation: Stability and Super-resolution, 2014. 
L= round(fc/2);
MusicHankel = @(y)hankel(y(1:L),y(L:fc));
[U, S, V] = svd(MusicHankel(y),0);

% ind=find(S>tol);
% N=ind(end)+1;
U_bot = U(:,N+1:end);

B = [];
for j=1:size(U_bot,2)
    u = U_bot(:,j);
    v = conj(u)';
    B(:,j)=conv(u,v);
end
C=sum(B,2);


R_n = roots(C(end:-1:1));
% figure
% plot(exp(2i*pi*(0:0.01:1)), 'k');
% hold on
% plot(R_n,'.')
R1_n = R_n(abs(R_n)<=1.000001); %find all the points on and inside the unit circle

[xi_es,I]=sort(mod(angle(R1_n),2*pi)/(2*pi)); %estimate possible positions
R1_n = R1_n(I);
% Keep only the best N ones. (closest to unit circle)
[~,I] = sort(abs(abs(R1_n)-1));
R1_n = R1_n(I);
xi_es = xi_es(I);

%check for double
ZZ = [angle(R1_n) abs(R1_n)];
[~,I]=sort(ZZ(:,1));
ZZ = ZZ(I,:);
ind_rem = find (abs(ZZ(1:end-1,1)-ZZ(2:end,1))<1e-3 &  abs(ZZ(1:end-1,2)-ZZ(2:end,2))<1e-3 );
xi_es(I(ind_rem))=[];
R1_n(I(ind_rem))=[];


xi_es =xi_es(1:N);
%plot(R1_n((1:N)),'.r')
Phi=@(x) exp(-2i*pi*(0:fc)'*x);
c_es = real(Phi(xi_es')\y);
%lsqnonneg(A,ob)