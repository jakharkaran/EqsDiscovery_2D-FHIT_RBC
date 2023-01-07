function uF = spectralFilter_same_size_x(u,NNX)
%% Sharp Spectral Filter for 2D Rectangular grids applied in x domain
% (implemented in spectral domain)

% The code finds and keeps the data from wave numbers below the cut-off
% wavenumber (smaller wave numbers)
% The data of higher wavenumbers is equated to zero (above cut-off
% wavenumber)

%% Input
% u: input flow filed (usually DNS): rectangular matrix
% NNX: 1x1 array

%% Output 
% uF: filtered flow field: square matrix of size(u)

u_hat = fft(u,[],2);
N = size(u,2);
filterSize = NNX/2;
u_hat(:,filterSize+2:N-filterSize+1)=0;
uF = real(ifft(u_hat,[],2));
end

