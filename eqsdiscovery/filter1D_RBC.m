function [Uf_c] = filter1D_RBC(U_DNS,filterType,coarseGrainingType, Delta, N_LES)
%% filter1D Filters and coarse grains 1D Rectangular grids (in x domain) 
% Use Case: Tuned for Rayleigh-Bernard Convection (RBC)

% filterType: gaussian
% 'gaussian': Gaussian Filter (Implemented in spectral domain)
% 'boxSpectral', 'boxPhysical', 'box': Box or Top hat filter 
%  (Implemented both in spectral and physical domain
% 'spectral' : Sharp Spectral filter (Implemented in spectral domain)

% Coarse graining
% 'spectral' : Removing high wave-numbers in spectral domain
% 'physical' : subsampling in physical domain

N_DNS = flip(size(U_DNS));

%% DNS grid RBC
N = 399;
NY = N + 1;
Lx = 6*pi;
kx = (2*pi/Lx)*[0:(N_DNS(1)/2) (-N_DNS(1)/2+1):-1];
Kx = ones(NY,1)*kx;

if N_DNS == N_LES
    coarseGrainingType = 'none';
end

%% Filtering

if strcmp(filterType,'gaussian')

    Gk = exp(-Kx.^2*Delta^2/24);
    Uf_hat = Gk.*fft(U_DNS,[],2);
    Uf = (ifft(Uf_hat,[],2));

elseif strcmp(filterType,'box') || strcmp(filterType,'boxSpectral')

    Gk = sin(0.5 * Kx* Delta)./(0.5*Kx*Delta);
    Gk(:,1) = 1;
    Uf_hat = Gk.*fft(U_DNS,[],2);
    Uf = (ifft(Uf_hat,[],2));

elseif strcmp(filterType,'gaussian+box') || strcmp(filterType,'box+gaussian') || ...
        strcmp(filterType,'gaussian+boxSpectral') || strcmp(filterType,'boxSpectral+gaussian')

    Gkbox = sin(0.5 * Kx* Delta)./(0.5*Kx*Delta);
    Gkbox(:,1) = 1;
    Gkgaussian = exp(-Kx.^2*Delta^2/24);
    Gk = Gkgaussian .* Gkbox;

    Uf_hat = Gk.*fft(U_DNS,[],2);
    Uf = (ifft(Uf_hat,[],2));

elseif strcmp(filterType,'spectral')

    Uf = spectralFilter_same_size_x(U_DNS,N_LES(1));

elseif strcmp(filterType,'none')
    Uf = U_DNS;
end


%% Coarse Graining
if strcmp(coarseGrainingType, 'spectral')

    Uf_c = spectralCoarse_x(Uf,N_LES(1));

elseif strcmp(coarseGrainingType, 'physical')
% Subsampling in physical space
    error('Code for subsampling in physical space is missing')

elseif strcmp(coarseGrainingType, 'none')
    Uf_c = Uf;
end

end