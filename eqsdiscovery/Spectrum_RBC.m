function [A_hat_mean] = Spectrum_RBC(A)

% Angled averaged spectrum for RBC
% Input: 2D Rectangular matrix (400 x N_LES)
% Output: 1D array

%% Checking the validity of input data
sizeA = size(A);

if length(sizeA) ~= 2
    error('Input flow variable is not 2D. Please in input 2D matrix.');
end

N_LES = sizeA(2);

% meanA = mean(A);
% meanA = ones(size(A,1),1)*meanA;
% Aprime = A - meanA;
Aprime = A;

% fft
A_hat = real(fft(Aprime,[],2));

% Averaging spectra over vertical direction
A_hat_mean = mean(A_hat(:,1:N_LES/2),1);

end
