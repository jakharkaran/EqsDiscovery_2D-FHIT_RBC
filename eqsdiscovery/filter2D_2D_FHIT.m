function [Uf_c] = filter2D_2D_FHIT(U_DNS,filterType,coarseGrainingType, Delta, N_LES)
% filter2D Filters and coarse grains 2D square grids
% Use Case: 2D Forced Homogenous Isotropic Tubulence (2D-FHIT)

%% Input

% U_DNS: Input DNS flow field: Square Matrix

% filterType: Characters
% 'gaussian': Gaussian Filter (Implemented in spectral domain)
% 'boxSpectral', 'box': Box or Top hat filter implemented in spectral domain
% 'spectral' : Sharp Spectral filter (Implemented in spectral domain)
% 'gaussian+boxSpectral', 'boxSpectral+gaussian', 'gaussian+box', 'box+gaussian': 
%  Gaussian and Box (spectral) Filters

% coarseGrainingType: Characters
% 'spectral' : Removing high wave-numbers in spectral domain
% 'physical' : subsampling in physical domain

% Delta: filter width: positive real number

% N_LES: size of coarse (LES) grid: 2x1 array

%% Output
% Uf_c: filtered and coarse grained flow filed: square matric of size N_LES

N_DNS = size(U_DNS);

%% DNS grid 2D Turbulence

Lx = 2*pi; % Length of domain
dx = Lx/N_DNS(1);
% x = linspace(0,Lx-dx,N_DNS(1));
kx = (2*pi/Lx)*[0:(N_DNS(1)/2) (-N_DNS(1)/2+1):-1];

% Making meshgrid
[Ky,Kx] = meshgrid(kx,kx);
Ksq = Kx.^2 + Ky.^2;

% Coarse grid
% kkx = (2*pi/Lx)*[0:(N_LES(1)/2) (-N_LES(1)/2+1):-1];
% [KKy,KKx] = meshgrid(kkx,kkx);

% No coarse graining
if N_DNS == N_LES
    coarseGrainingType = 'none';
end

%% Filtering

if strcmp(filterType,'gaussian')

    Gk = exp(-Ksq*Delta^2/24); % transfer function
    Uf = ifft2(Gk.*fft2(U_DNS));

elseif strcmp(filterType,'box') || strcmp(filterType,'boxSpectral')

    Gkx = sin(0.5 .* Kx * Delta)./(0.5.*Kx*Delta);
    Gky = sin(0.5 .* Ky * Delta)./(0.5.*Ky*Delta);
    Gkx(1,:) = 1;
    Gky(:,1) = 1;
    Gk = Gkx .* Gky; % transfer function

    Uf = ifft2(Gk.*fft2(U_DNS));

elseif strcmp(filterType,'gaussian+box') || strcmp(filterType,'box+gaussian') || ...
        strcmp(filterType,'gaussian+boxSpectral') || strcmp(filterType,'boxSpectral+gaussian')
    Gkx = sin(0.5 .* Kx * Delta)./(0.5.*Kx*Delta);
    Gky = sin(0.5 .* Ky * Delta)./(0.5.*Ky*Delta);
    Gkx(1,:) = 1;
    Gky(:,1) = 1;
    Gkbox = Gkx .* Gky; % box filter transfer function

    Gkgaussian = exp(-Ksq*Delta^2/24); % gaussian filter transfer function
    Gk = Gkgaussian .* Gkbox; % transfer function
    Uf = ifft2(Gk.*fft2(U_DNS));

elseif strcmp(filterType,'spectral_circle') || strcmp(filterType,'spectral')
    Uf = spectralFilter_circle_same_size_2D(fft2(U_DNS),N_LES(1));

elseif strcmp(filterType,'none')
    Uf = U_DNS;
end


%% Coarse Graining

if strcmp(coarseGrainingType, 'spectral')
% Coarse graining in spectral space
    Uf_c = spectralCoarse_square_2D(fft2(Uf),N_LES(1));

elseif strcmp(coarseGrainingType, 'physical')
% Subsampling (coarse graining) in physical space

    xgrid = (1:N_LES(1))*(N_DNS(1)/N_LES(1)) - (N_DNS(1)/N_LES(1))/2;
    ygrid = (1:N_LES(2))*(N_DNS(2)/N_LES(2)) - (N_DNS(2)/N_LES(2))/2;

    Uf_c = Uf(xgrid,ygrid);

elseif strcmp(coarseGrainingType, 'none')
    Uf_c = Uf;
end

end