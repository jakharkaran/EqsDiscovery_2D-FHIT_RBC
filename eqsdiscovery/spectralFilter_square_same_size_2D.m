function uF = spectralFilter_square_same_size_2D(u,NNX)
%% Sharp Spectral Filter for 2D square grids 
% (implemented in spectral domain)

% The code finds and keeps the data from wave numbers below the cut-off
% wavenumber (smaller wave numbers)
% The data of higher wavenumbers is equated to zero (above cut-off
% wavenumber)

% This is implemented in square domain ({Kx,Ky > Kc} = 0)
% Kx, Ky : wave number in x and y domain respectively
% Kc : Cut-off wave number

%% Input
% u: input flow filed (usually DNS): square matrix
% NNX: 1x1 array

%% Output 
% uF: filtered flow field: square matrix of size(u)

N = size(u,1);
filterSize = NNX/2;
u(filterSize+2:N-filterSize+1,:)=0;
u(:,filterSize+2:N-filterSize+1)=0;
uF = real(ifft2(u));
% uF = u;
end