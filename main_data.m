%%
%-------------------------------------------------------------------------%
%                             Ph.D. THESIS                                %
%           LINEAR UNMIXING IN HYPERSPECTRAL IMAGERY USING ADMM           %
%                             VARIABILITY                                 %
%-------------------------------------------------------------------------%
%% File data
% File : main_data.m
% Author : P.A. Thouvenin (22/10/2014)
% Last modified : 05/11/2014
clc, clear all, close all;
%-------------------------------------------------------------------------%
%% Synthetic data parameters
K        = 3;    % Number of endmembers in the scene
H        = 128;  % Height of the synthetic image
W        = 64;   % Width    ------------------
N        = H*W;  % Number of pixel
SNR      = 15;   % SNR of the data (dB)
sigma2   = 300; % Smoothness parameter for the synthetic data
cutoff   = 0.7;  % Maximum abundance value (no pure pixels)
position = [20 5; 110 60];
coeffvar = [0.1,0.1,0.25,0.25]; % variability coefficient applied to each of the four subpart of the image (image divided into four tiles)
H1 = H/2;
W1 = W/2;

%-------------------------------------------------------------------------%
%% Synthetic data generation
disp('Generating synthetic data...');
[data,M_var,M,dM,A] = generate_data(H,H1,W,W1,K,SNR,coeffvar,sigma2,position,cutoff);
disp('... DONE');
disp('---------------------------------------------------------------------------');

% Data size
[H,W,L] = size(data);

% Theoretical abundance map and spectra
figure('Name','Theoretical endmembers','NumberTitle','Off');
plot(M)

A1 = permute(reshape(A',W,H,K),[2 1 3]); 
for k = 1:K
    % Results display
    figure('Name',['Theoretical abundance and spectra - endmember ' num2str(k)],'NumberTitle','Off');
    % Abundance map
    subplot(2,1,1);
    imagesc(A1(:,:,k));
    colormap('hot');
    colorbar;
    title('Abundance map');
    % Endmember
    subplot(2,1,2);
    plot(M(:,k)) %,'Color',color(k,:)
    title('Endmembers');
end
clear A1;
%-------------------------------------------------------------------------%