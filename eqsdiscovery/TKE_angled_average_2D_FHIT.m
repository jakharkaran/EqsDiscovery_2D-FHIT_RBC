function [TKE, E] = TKE_angled_average_2D_FHIT(U,V)

% Angled averaged energy spectrum for 2D_FHIT

% Input: 2D velocity inputs (U, V): Square Matrices

% Output
% TKE: Angled average Turbulent Kinetic Energy (TKE): 1D array
% E: Kinetic Energy in physical domain

%% Checking the validity of input data

sizeU = size(U);
sizeV = size(V);
if length(sizeU) ~= 2 || length(sizeV) ~= 2
    error('Input velocity is not 2D. Please in put 2D velocity.');
end
if sizeU(1) ~= sizeU(2) || sizeV(1) ~= sizeV(2)
    error('Dimesion mismatch for velocity. Both dimensions of velocity should be equal.');
end
if length(U) ~= length(V)
    error('Velocity do not have equal dimensions.');
end

NX = length(U);
Lx = 2*pi;
kx = (2*pi/Lx)*[0:(NX/2) (-NX/2+1):-1];

[Ky,Kx] = meshgrid(kx,kx);
Ksq = Kx.^2 + Ky.^2;
invKsq = 1./Ksq;
invKsq(1,1) = 0;
Kabs = sqrt(Ksq);

U = U - mean(U,'all');
V = V - mean(V,'all');


E = 0.5*(U.^2 + V.^2);
Ehat = 0.5*(fft2(U.^2 + V.^2));

% Angle averaged energy spectrum
%% Why are you using this particular length of an array?

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
TKE = tke/(NX^2);

