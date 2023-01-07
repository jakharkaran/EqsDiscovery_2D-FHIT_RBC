function [Tdash, TnameOut] = derivative_2D_FHIT(T,order,Tname)
%% Calculate spatial derivatives for 2D_FHIT
% Derivatives are calculated in spectral space
% Boundary conditions are periodic in x and y spatial dimensions
% Length of domain: 2*pi

% [Tdash, 'Txxy'] = derivative_2D_FHIT(T, 2, 1, 'T')

%% Input
% T: Input flow field: Square Matrix NxN
% order [orderX, orderY]: Respective order of derivatives in x and y spatial dimensions: >=0 Integers
% Tname: Name of the input flow field: Characters 'T', 'Psi', 'tau', 'U',...

%% Output
% Tdash: derivative of the flow field T: Square Matrix NxN
% TnameOut: Name of the output derivative of flow field: Characters

orderX = order(1);
orderY = order(2);

% Validating user input
if orderX < 0 || orderY < 0
    error('Order of derivatives must be 0 or positive');
elseif orderX == 0 && orderY == 0
    error('Both order of derivatives are 0, atleast one of them should be positive');
end

% orderX 0,1,2,3,4...
% orderY,0,1,2,3,4...

Ngrid = size(T);

Lx = 2*pi; % Length of domain
kx = (2*pi/Lx)*[0:(Ngrid(1)/2) (-Ngrid(1)/2+1):-1];

% Making meshgrid
[Ky,Kx] = meshgrid(kx,kx);

% Calculating derivatives in spectral space

T_hat = fft2(T);
Tdash_hat = ((1i*Kx).^orderX) .* ((1i*Ky).^orderY) .*T_hat;
Tdash = real(ifft2(Tdash_hat)) ;

% Naming the variable
TnameOut = Tname;
if orderX > 0
    for count = 1:orderX
        TnameOut = [TnameOut 'x'];
    end
end
if orderY > 0
    for count = 1:orderY
        TnameOut = [TnameOut 'y'];
    end
end

end