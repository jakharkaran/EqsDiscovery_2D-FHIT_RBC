function [A_hat_mean] = TKE_RBC(U,V)

% TKE spectrum for RBC
% Input: 2D Rectangular matrix (x x y) spatial dimension eg: N_LES (N_DNS) x 400
% Output: 1D array

%% Checking the validity of input data
sizeU = size(U);
sizeV = size(V);

if length(sizeU) ~= 2 && length(sizeV) ~= 2
    error('Input flow variable is not 2D. Please in input 2D matrix.');
end

N_LES = sizeU(1);

meanU = mean(U,'all');
meanV = mean(V,'all');

Uprime = U - meanU;
Vprime = V - meanV;

% Energy
E = Uprime.^2 + Vprime.^2;

% fft
E_hat = fft(E,[],1);

% Averaging spectra over vertical direction
A_hat_mean = mean(abs(E_hat(1:N_LES/2,:)),2)/N_LES;

end
