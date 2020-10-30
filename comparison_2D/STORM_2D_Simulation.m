function [im, emitterList,XX,YY] = STORM_2D_Simulation(density, pixelsize, ...
    width, height, g_noise, photons, sigma)
% PURPOSE:
% Generate simulated single snapshot 2D STORM image.
%---------------------------------------------------
% USAGE:
% [im, emitterList] = STORM_2D_Simulation(density, pixelsize, ...
%	width, height, g_noise, photons, sigma)
%---------------------------------------------------
% INPUTS:
% density:      emitter density (fluorophores per um^2)
% pixelsize:    camera pixel size (um)
% width:        image width (pixel)
% height:       image height (pixel)
% g_noise:      gaussian noise level (photons/pixel)
% photons:      average photon number
% sigma:        standard deviation of system PSF (pixel)
%---------------------------------------------------
% OUTPUTS:
% im:           simulated image (noised)
% emitterList:  emitter coordinates
%---------------------------------------------------

myseed=RandStream.create('mrg32k3a', 'NumStreams', 2^10, 'seed', 'shuffle');
RandStream.setGlobalStream(myseed);

xx = linspace(1, width, width)';
yy = linspace(1, height, height)';
[XX, YY] = meshgrid(xx, yy);

%% Generate the emitterlist
margin = 3;
% average active fluorophores per frame
P = density * (width - 2 * margin - 1) ^ 2 * pixelsize ^ 2;
% round P to integer randomly
if rand() < (P - floor(P))
    P = floor(P) + 1;
else
    P = floor(P);
end

% generate emitters location randomly
emitterList = rand(P, 2);
% emitterList(:, 1) = [6.7 6.2];
% emitterList(:, 2) = [6.7 5];

emitterList(:, 1)=emitterList(:, 1) * (width - 2 * margin - 1) + margin + 1;
emitterList(:, 2)=emitterList(:, 2) * (height - 2 * margin - 1) + margin + 1;
%emitterList = emitterList';
%% Generate the image
I = zeros(height, width);
% iterate through all emitters
for j = 1 : size(emitterList, 1)
    r0_vec = emitterList(j, :);
    I = I + photons / 4 *(erf((XX + 0.5 - r0_vec(1)) / sqrt(2) / sigma) - erf((XX - 0.5 - r0_vec(1)) / sqrt(2) / sigma)) .*(erf((YY + 0.5 - r0_vec(2)) / sqrt(2) / sigma) -erf((YY - 0.5 - r0_vec(2)) / sqrt(2) / sigma));
end
% Poisson noise
 I_noise = uint16(I);
 I_noise = imnoise(I_noise, 'poisson');
% % Gaussian noise
 I_noise = double(I_noise) + g_noise .* randn(width, height);
 I_noise = uint16(I_noise);
%im = double(I_noise);
im = I;
%plot(im-im1)