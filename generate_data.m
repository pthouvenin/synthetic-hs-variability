function [data,M_var,M,dM,A] = generate_data(H,H1,W,W1,K,SNR,coeff_var,sigma2,position,cutoff)
% DATA_VAR_SMOOTH_SYNTH generates a smooth synthetic hyperspectral 
% datacube from a given number of variable endmembers.
% > type         Abundance spatial configuration ('smooth'|'nsmooth)
% > H    (1,1)   Height of the datacube
% > W    (1,1)   Width of the datacube
% > P    (1,1)   Number of endmember classes in the scene
% > SNR  (1,1)   SNR (dB) of the data generated (gaussian noise assumption)
% > coeffvar     Variability of each class from the initial endmember (1|4)
% > sigma2       Spatial smoothness parameter
% > position     Endmember position in the image (max. abundance)
% > cutoff       Maximal abundance value
%
% < data  (H,W,L) Generated synthetic hyperspectral datacube
% < M_var (1,N)   Perturbed endmembers (cell)
% <  M    (L,P)   Theoretical endmembers
% < dM    (1,N)   Theoretical perturbation (cell)
% <  A    (P,H*W) Abundance matrix
%% Initial test
if (H1 >= H) || (W1 >= W)
    error('Invalid H1 or W1 : H > H1 and W > W1 must be repsected')
end

%% Initialisation
load spectres.txt
spectra = spectres(:,2:end);
clear spectres

Nsp = size(spectra,2);
M = spectra(1:2:end,randperm(Nsp,K));
L = size(M,1);
N = H*W;
dM = cell(1,N);

%% Abundance generation
A = smooth_abundance(sigma2,H,W,K,cutoff,position);

%% Variability and data generation
M_var = cell(1,N);
Mvar  = zeros(L,K);
data  = zeros(L,N);

for h = 0:H-1
    for w = 1:W
        n = w + h*W;
        if (h <= H1-1)
            if (w <= W1)
                coeffvar = coeff_var(1);
            else
                coeffvar = coeff_var(2);
            end
        else
            if (w <= W1)
                coeffvar = coeff_var(3);
            else
                coeffvar = coeff_var(4);
            end
        end
        
        if ~isequal(coeffvar,0)
            for k = 1:K
                Mvar(:,k) = generate_variability(M(:,k),L,coeffvar);
            end
            M_var{n} = Mvar;
            dM{n} = M - Mvar;   % computation of the perturbation
            data(:,n) = Mvar*A(:,n);        
        else
            M_var{n} = M;
            dM{n} = zeros(L,K);
            data(:,n) = M*A(:,n);
        end
    end
end

data = permute(reshape(data',W,H,L),[2 1 3]);

%% Noising
Py = mean(data(:).^2);
sigma2 = Py/10^(SNR/10);
data = abs(data + randn(H,W,L)*sqrt(sigma2));

end
