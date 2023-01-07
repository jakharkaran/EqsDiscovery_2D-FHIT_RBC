function [Tdash, TnameOut] = derivative_RBC(T,order,Tname)
% Calculate spacial derivatives for RBC
% T: Input flow field
% orderX, orderY: Respective order of derivatives in x and y spacial dimensions
% Tname: Character input naming the derivative

orderX = order(1);
orderY = order(2);

if orderX < 0 || orderY < 0
    error('Order of derivatives must be 0 or positive');
elseif orderX == 0 && orderY == 0
    error('Both order of derivatives are 0, atleast one of them should be positive');
end

% orderX 0,1,2,3,4
% orderY,0,1,2,3,4

sizeT = size(T);

NX = sizeT(2); % Grid size in x
N = 399;
NY = N + 1; % Grid size in y

% Derivatives

T_hat = fft(T,[],2);

% chebishev grid in y spatial domain
if orderY > 0

    Ty_hat = zeros([NY,NX]);
    Dy = -cheb(N); % First derivative in y

    Ty_hat(:,1) = (Dy)^orderY*T_hat(:,1);

    for m = 2:NX/2       
        Ty_hat(:,m) = (Dy)^orderY*T_hat(:,m);
        Ty_hat(:,NX-(m-2)) = conj(Ty_hat(:,m));
    end
    Ty_hat(:,NX/2+1) = (Dy)^orderY*T_hat(:,NX/2+1);
    Ty = real(ifft(Ty_hat,[],2));
end

% Peridic in x spatial domain
if orderX > 0

    Lx = 6*pi;
    kx = (2*pi/Lx)*[0:(NX/2) (-NX/2+1):-1];

    Tx_hat = zeros([NY,NX]);
    Tx_hat(:,1) = 0;

    if orderY > 0
        T_hat = Ty_hat;
    end

    for m = 2:NX/2   
        Tx_hat(:,m) =(1i*kx(m)).^orderX*T_hat(:,m);
        Tx_hat(:,NX-(m-2)) = conj(Tx_hat(:,m));
    end
    Tx_hat(:,NX/2+1) = (1i*kx(NX/2+1)).^orderX*T_hat(:,NX/2+1);
    Tx = real(ifft(Tx_hat,[],2));
end

if orderX > 0
    Tdash = Tx;
else
    Tdash = Ty;
end

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