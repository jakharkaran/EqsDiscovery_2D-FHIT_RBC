function [A_hat_mean] = spectrum_RBC(A)

% Spectrum for RBC
% Input: 2D Rectangular matrix (x x y) spatial dimensions eg: N_LES (N_DNS) x 400
% Output: 1D array

%% Checking the validity of input data
sizeA = size(A);

if length(sizeA) ~= 2
    error('Input flow variable is not 2D. Please in input 2D matrix.');
end

N_LES = sizeA(1);

meanA = mean(A,'all');
Aprime = A - meanA;

% fft
A_hat = fft(Aprime,[],1);

% Averaging spectra over vertical direction
A_hat_mean = mean(abs(A_hat(1:N_LES/2,:)),2)/N_LES;

end
