%This code has been modified by A. Koulouri in order to compare the standard STORM (https://users.ece.cmu.edu/~yuejiec/publications.html publication:)  with
%the adaptive scheme

close all; clear
addpath(genpath(cd)) 
%% Set up the parameters
density = 10;           % emitter density (fluorophores per um^2)
pixelsize = 0.075;      % camera pixel size (um)
width = 18;             % image width (pixel)
height = 18;            % image height (pixel)
g_noise = 0;            % gaussian noise level (photons/pixel) We add
photons = 500;         % average photon number
sigma = 1.5;            % standard deviation of system PSF (pixel)

%% Generate measurement matrix
disp('Generating measurement matrix ...');
div1 = 5; % zoom factor
% measurement matrix for CS
A1 = STORM_2D_Gen_Meas_Mat(width, height, div1, sigma);

% measurement matrix for ADCG
div2 = 4;
A2 = STORM_2D_Gen_Meas_Mat(width, height, div2, sigma);

%% Generate simulated image
disp('Generating simulated image ...');
[im, emitterList,xm,ym] = STORM_2D_Simulation(density, pixelsize, ...
    width, height, g_noise, photons, sigma);
SNR = 60; %Noise level
n = randn(size(im));
scale = (norm(im(:))^2/norm(n(:))^2)* 10^(-SNR/10);
im = im + sqrt(scale)*n;
im(im<1e-10) = 0;

%% Compressed sensing on the image (SC-STORIM)
thresh = 1;
ns = 1.5;
di_wid = 1;
t1 = tic;
[clustCent1, zm] = STORM_2D_CS(A1, im, height, width, div1, thresh, ns, di_wid);
disp(['CS processing time: ', num2str(toc(t1)), ' seconds']);


 ns = 0.9;
 thresh = 10;
 t1 = tic;
 [clustCent2, I_list2] = TVSTORM_2D(A2, im, height, width, sigma, div2, ns, thresh);
 disp(['TVSTORM processing time: ', num2str(toc(t1)), ' seconds']);

%% Adaptive compressed sensing
Lim = 0.15;
size_comp_grid = round(width/3); %select initial grid
QuietPlot =1;
[xi_es,c_es]= AS_alg(im(:),xm,ym,sigma,1,width,height,size_comp_grid,Lim,QuietPlot,emitterList(:,1),emitterList(:,2));



%%Display
figure;
im = reshape(im, height, width);
imagesc(im); colormap gray; axis equal; axis off; hold on;
plot(emitterList(:, 1), emitterList(:, 2), 'kx', 'MarkerSize', 7, 'LineWidth', 1.5);
plot(clustCent1(:, 1), clustCent1(:, 2), 'y*', 'MarkerSize', 7, 'LineWidth', 1.5);
plot(clustCent2(:, 1), clustCent2(:, 2), 'bo', 'MarkerSize', 7, 'LineWidth', 1.5);
%plot(emitterList(:, 1), emitterList(:, 2), 'kx', 'MarkerSize', 7, 'LineWidth', 1.5);
plot( xi_es(:,1), xi_es(:,2), 'ro', 'MarkerSize', 7, 'LineWidth', 1.5);

legend('true','CSSTORM','TVSTORM','Adaptive L1')


axis off

% psize = get(gcf,'PaperSize');
% wd = 7;
% hg =7; 
% lf = (psize(1)-wd)/2;bt = (psize(2)-hg)/2;
% set(gcf,'PaperPosition',[lf bt wd hg]);
%  print('-dpng','-r600',['test_4'])
%  print('-depsc2','-r700',['test_4'])
