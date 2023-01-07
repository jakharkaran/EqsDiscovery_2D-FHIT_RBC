function uF = spectralFilter_circle_2D(u,NNX)
%% Coarse graining followed by spectral filter (both in spectral domain)

% The code finds and keeps the data from wave numbers below the cut-off
% wavenumber (smaller wave numbers) and stores them in lower dimensional
% matrix (coarse graining). This is followed by applying a spectral filter
% in coarse grained grid.

% For spectral filter applied in spectral domain:
% filtering followed by coarse graining or vice versa will provide with same result

%% Input
% u: input flow filed (usually DNS): square matrix
% NNX: side length of output flow field: 1x1 array

%% Output 
% uF: coarse grained flow field: square matrix of size: NNX x NNX

N = size(u,1);
filterSize = NNX/2;
u = fftshift(u);
u_hat = fftshift(u((N/2)-filterSize+2:(N/2)+1+filterSize,(N/2)-filterSize+2:(N/2)+1+filterSize))/((N/2)/filterSize)^2;
u_hat = circshift(u_hat,1,1);
u_hat = circshift(u_hat,1,2);

u_hat(filterSize+1,:)=0;
u_hat(:,filterSize+1)=0;


indX = zeros(NNX,NNX);
indY = zeros(NNX,NNX);

sizeIndX = size(indX);
sizeIndY = size(indY);

[indY, indX] = meshgrid(1:sizeIndX(1),1:sizeIndY(1)); 
    
dist = zeros(NNX,NNX);

subTractX = 1;
subTractY = 1;
dist(1:filterSize,1:filterSize) =  sqrt((subTractX - indX(1:filterSize,1:filterSize)).^2 ...
    + (subTractY - indY(1:filterSize,1:filterSize)).^2);

subTractX = 1;
subTractY = NNX;
dist(1:filterSize,NNX-filterSize+1:NNX) =  sqrt((subTractX - indX(1:filterSize,NNX-filterSize+1:NNX)).^2 ...
    + (subTractY - indY(1:filterSize,NNX-filterSize+1:NNX)).^2);

subTractX = NNX;
subTractY = 1;
dist(NNX-filterSize+1:NNX,1:filterSize) =  sqrt((subTractX - indX(NNX-filterSize+1:NNX,1:filterSize)).^2 ...
    + (subTractY - indY(NNX-filterSize+1:NNX,1:filterSize)).^2);

subTractX = NNX;
subTractY = NNX;
dist(NNX-filterSize+1:NNX,NNX-filterSize+1:NNX) =  sqrt((subTractX - indX(NNX-filterSize+1:NNX,NNX-filterSize+1:NNX)).^2 ...
    + (subTractY - indY(NNX-filterSize+1:NNX,NNX-filterSize+1:NNX)).^2);

dist(dist>filterSize) = 0;
dist(dist>0) = 1;
dist(1,1) = 1;
dist(1,NNX) = 1;
dist(NNX,1) = 1;
dist(NNX,NNX) = 1;

u_hat = u_hat.* dist;

uF = real(ifft2(u_hat));
% uF = u_hat;
end