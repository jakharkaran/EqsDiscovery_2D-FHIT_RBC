function [TKE] = spectrum_angled_average_2D_FHIT(U)

% Angled averaged spectrum for 2D flow variables
% Input: 2D Square matrix
% Output: 1D array

%% Checking the validity of input data
sizeU = size(U);

if length(sizeU) ~= 2
    error('Input flow variable is not 2D. Please in input 2D matrix.');
end
if sizeU(1) ~= sizeU(2)
    error('Dimension mismatch for flow variable. Flow variable should be a square matrix.');
end

NX = length(U);
Lx = 2*pi; % Length of domain
kx = (2*pi/Lx)*[0:(NX/2) (-NX/2+1):-1];

[Ky,Kx] = meshgrid(kx,kx);
Ksq = Kx.^2 + Ky.^2;
invKsq = 1./Ksq;
invKsq(1,1) = 0;
Kabs = sqrt(Ksq);

% Removing mean
U = U - mean(U,'all');

% E = U;
Ehat = fft2(U);

% Angle averaged energy spectrum
arr_len = floor(0.5*sqrt((NX*NX + NX*NX)));
n = arr_len + 1;
kplot = linspace(2,n,n);

eplot = zeros(n,1);

% 
es = pi*(abs(Ehat))./Kabs;
% es = abs(wor_hat);

for k = 2:n 
    eplot(k) = 0;
    ic = 0;
    for i = 2:NX
        for j = 2:NX
            kk = sqrt(kx(i).^2 + kx(j).^2);
            if (kk >= (k - 0.5) && kk < (k+0.5))
                ic = ic+1;
                eplot(k) = eplot(k) + es(i,j);
            end
        end
    end
%     eplot(k) = eplot(k) / ic;
end

%%

tke = mean(eplot,2);
TKE = tke/NX^2;
