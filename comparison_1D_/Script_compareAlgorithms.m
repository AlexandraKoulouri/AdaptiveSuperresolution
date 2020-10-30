%this script was created by A. Koulouri 19.10.2020
%We run a test with 5 peaks using Dirichlet kernel
%we compare the deconvolution results using the adaptive grid l1-norm,
%MUSIC (Single-Snapshot Spectral Estimation) and SDP(Super-resolution from noisy data" 
% by E. Candes and C. Fernandez-Granda)


clear
close all

close all; clear
addpath(genpath(cd)) 

fc = 18; %number of frequencies
Q=1024*4;
z = (0:Q)'/Q; % measurement/observation grid size

%locations of the peaks and amplitudes
xi =[0.13 0.453 0.15 0.7 0.9];
N=length(xi);
c =[0.7 -0.8 0.9 1 0.9];


Phi=@(x) exp(-2i*pi*(0:fc)'*x);

sigma = 1e-1;
nn =sigma* randn(fc+1,1);
w = (randn(fc+1,1)+1i*randn(fc+1,1)); w = w/norm(w);
y0 = Phi(xi)*c(:);
y_fourier = y0+nn; %sigma*norm(y0)*w; Signal in the Fourier Domain

%in the time domain the signal(take the inverse Fourier transform)
PhiT =@(x)exp(2i*pi*x*(0:fc))/(fc+1);
y = real(PhiT(z)*y_fourier); %Signal in Time domain


[xi_es_L1,c_es_L1] = AS_1D(y,30,z,fc,xi,c,0.5/fc,0.1); %solution using Adaptive L1-norm in 1D

[xi_es_music,~] = RootMUSIC(y_fourier,5,fc); %solution using MUSIC
A = DirichletKernel(fc,xi_es_music',z);
c_es_music =  l1_ls(A,y,0.1*max(A'*y)); %amplitude based on estimated locations from MUSIC

lambda= 1*sigma*sqrt(fc*log(fc));
[xi_es_SDP,c_es_SDP]=SDP_solver([conj(y_fourier(end:-1:2));(y_fourier)],fc,y,z,lambda);%solution using 



figure;
stem(xi,c,'rx','DisplayName','true peak');
hold on
plot(xi_es_L1, c_es_L1,'ok','markersize',5,'color',[0 0 0],'DisplayName','adaptive L1')
plot(xi_es_music, c_es_music,'g*','markersize',5,'DisplayName','MUSIC')
plot(xi_es_SDP, c_es_SDP,'b+','markersize',5,'DisplayName','SDP')
%stem(xi,c,'rx','linewidth',1,'markersize',5.5)
plot(z,y,'color',[0.1 0.5 0.7],'linewidth',0.6,'DisplayName','observation')
 legend('true','Adaptive L1','MUSIC','SDP')