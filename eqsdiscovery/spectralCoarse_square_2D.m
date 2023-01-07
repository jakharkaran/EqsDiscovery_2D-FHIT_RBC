function uF = spectralCoarse_square_2D(u,NNX)
%% Coarse graining in spectral domain

% The code finds and keeps the data from wave numbers below the cut-off
% wavenumber (smaller wave numbers) and stores them in lower dimensional
% matrix (coarse graining)

% The data is implemented in square domain ({Kx,Ky > Kc} = 0)
% Kx, Ky : wave number in x and y domain respectively
% Kc : Cut-off wave number

%% Input
% u: fft2() of input flow filed (usually DNS): square matrix
% NNX: side length of output flow field: 1x1 array

%% Output 
% uF: coarse grained flow field: square matrix of size: NNX x NNX

N = size(u,1)/2;
filterSize = NNX/2;
u = fftshift(u);
u_hat = fftshift(u(N-filterSize+2:N+1+filterSize,N-filterSize+2:N+1+filterSize))/(N/filterSize)^2;
u_hat = circshift(u_hat,1,1);
u_hat = circshift(u_hat,1,2);

u_hat(filterSize+1,:)=0;
u_hat(:,filterSize+1)=0;
uF = real(ifft2(u_hat));
% uF = u_hat;
end
