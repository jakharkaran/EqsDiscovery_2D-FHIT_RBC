function uF = spectralFilter_circle_same_size_2D(u,NNX)
%% Sharp Spectral Filter for 2D square grids 
% (implemented in spectral domain)

% The code finds and keeps the data from wave numbers below the cut-off
% wavenumber (smaller wave numbers)
% The data of higher wavenumbers is equated to zero (above cut-off
% wavenumber)

% This is implemented in circular domain ({sqrt(Kx^2 + Ky^2) < Kc} = 0)
% Kx, Ky : wave number in x and y domain respectively
% Kc : Cut-off wave number

%% Input
% u: input flow filed (usually DNS): square matrix
% NNX: 1x1 array

%% Output 
% uF: filtered flow field: square matrix of size(u)

N = size(u,1);
filterSize = NNX/2;

indX = zeros(N,N);
indY = zeros(N,N);

sizeIndX = size(indX);
sizeIndY = size(indY);

[indY, indX] = meshgrid(1:sizeIndX(1),1:sizeIndY(1)); 

% Keeping the circular domain above threshold wave number and making rest
% of the data = 0
    
dist = zeros(N,N);

subTractX = 1;
subTractY = 1;
dist(1:filterSize,1:filterSize) =  sqrt((subTractX - indX(1:filterSize,1:filterSize)).^2 ...
    + (subTractY - indY(1:filterSize,1:filterSize)).^2);

subTractX = 1;
subTractY = N;
dist(1:filterSize,N-filterSize+1:N) =  sqrt((subTractX - indX(1:filterSize,N-filterSize+1:N)).^2 ...
    + (subTractY - indY(1:filterSize,N-filterSize+1:N)).^2);

subTractX = N;
subTractY = 1;
dist(N-filterSize+1:N,1:filterSize) =  sqrt((subTractX - indX(N-filterSize+1:N,1:filterSize)).^2 ...
    + (subTractY - indY(N-filterSize+1:N,1:filterSize)).^2);

subTractX = N;
subTractY = N;
dist(N-filterSize+1:N,N-filterSize+1:N) =  sqrt((subTractX - indX(N-filterSize+1:N,N-filterSize+1:N)).^2 ...
    + (subTractY - indY(N-filterSize+1:N,N-filterSize+1:N)).^2);

dist(dist>filterSize) = 0;
dist(dist>0) = 1;
dist(1,1) = 1;
dist(1,N) = 1;
dist(N,1) = 1;
dist(N,N) = 1;

u = u.* dist;
uF = real(ifft2(u));
% uF = u;
end