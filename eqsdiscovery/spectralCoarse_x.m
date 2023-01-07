function uF = spectralFilter_x(u,NNX)
%% Coarse graining in spectral domain 2D Rectangula grids applied in x domain
% Use case: Rayleigh-Bernard Convection

% The code finds and keeps the data from wave numbers below the cut-off
% wavenumber (smaller wave numbers) and stores them in lower dimensional
% matrix (coarse graining)

% The data is implemented in square domain ({Kx,Ky > Kc} = 0)
% Kx, Ky : wave number in x and y domain respectively
% Kc : Cut-off wave number

u = fft(u,[],2);
N = size(u,2)/2;
filterSize = NNX/2;
u = fftshift(u,2);
u_hat = fftshift(u(:,N-filterSize+2:N+1+filterSize),2)/(N/filterSize);
u_hat = circshift(u_hat,1,2);

u_hat(:,filterSize+1)=0;
uF = real(ifft(u_hat,[],2));
end