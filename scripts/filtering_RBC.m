addpath('../eqsdiscovery/');

%% Generate filtered and SGS flow fields for Rayleigh-Bernard Convection (RBC) on a rectangular grid

%% local
ICtype = 'train1';
dataType = 'NX2048NY400Ra4E7Pr7';
Ngrid = [128 256 512]; % Ngrid = [48 96 192 384];
Nfilter = Ngrid;
filterType = 'gaussian';
skipSnap = 1;
coarseGrainType = 'spectral';

% filterTypeArr = ["gaussian", "box", "gaussian+box"];
filterTypeArr = ["gaussian"];

%%
% ICType = 'train1'; % train1 test1

% dataType = 'NX2048NY400Ra1E6Pr100';
% NX2048NY400Ra1E6Pr100 v NX2048NY400Ra4E7Pr100 NX2048NY400Ra4E7Pr7

% Ngrid = [8,16,32,64,128, 256];                    % Nfilter = Ngrid for coarse grid
% Nfilter = 128;     % Nfilter = Ngrid for coarse grid (NNX != NX) 

% Select Filter type
% filterType = 'gaussian'; % gaussian box sharpSpectral

% Select filtered quantities to be saved | Equate to 1 to save 
saveVelocity = 1;   % Save filtered velocity (both U and V)
saveStreamfxnVorticity = 1;     % Save stream function and vorticity
saveTemperature = 1;
saveResidualStress = 1;     % Save filtered residual stress \bar(uv) - \bar(u)\bar(v)
saveResidualStressVor = 1;    % Save filtered residual stress of vorticity \bar(uVor) - \bar(u)\bar(Vor)
saveResidualStressT = 1; % Save filtered residual stress of Temperature (J) \bar(uT) - \bar(u)\bar(T)
savePIuv = 0;    % Save PIuv (Divergence of tau)
savePIvor = 1;    % Save PIvor (Curl of Divergence of tau or divergence of sigma)
savePIt = 1;     % Save PIt (Divergence of J)

% Conservative or non-conservative form of PI's can be calculated
PIconservative = 1; % Conservative form of PI will be calculated (default)
runCode = 'local'; % 'clusture': running the code on clusture, 'local': running the code locally

% Specific to DNS data
% skipSnap = 1;  % Save every skipSnap th filtered snapshot

switch runCode

    case 'clusture'
    
    if strcmp(dataType,'NX2048NY400Ra1E6Pr100') || strcmp(dataType,'NX2048NY400Ra4E7Pr100') ...
            || strcmp(dataType,'NX2048NY400Ra4E7Pr7')
        nosSnapshot = ones(1,17)*50; % number of snapshots in each file used to save DNS data
            if strcmp(ICtype,'test1')
            nosSnapshot = ones(1,2)*50;
            end
    end

    case 'local'

    if strcmp(dataType,'NX2048NY400Ra1E6Pr100') || strcmp(dataType,'NX2048NY400Ra4E7Pr100') ...
            || strcmp(dataType,'NX2048NY400Ra4E7Pr7')
        nosSnapshot = ones(1,1)*20; % number of snapshots in each file used to save DNS data
%             if strcmp(ICtype,'test1')
%             nosSnapshot = ones(1,2)*50;
%             end
    end
end

NX = 2048;
totalSnapshots = floor(sum(nosSnapshot)/skipSnap);  % total number of filtered snapshots saved

DATA_DIR = ['../data/RBC/' dataType '/DNS/' ICtype '/'];

dummy = 'dummy';

disp('******************************************************');
disp('---------------RBC----------------');
disp(['dataType                        = ' dataType]);
disp(['Filter Type                     = ' filterType]);
disp(['ICtype                          = ' ICtype ]);
disp(['N DNS                           = ' num2str(NX)]);
disp(['Ngrid                           = ' num2str(Ngrid)]);
disp(['Nfilter                         = ' num2str(Nfilter)]);
disp(['skipSnap                        = ' num2str(skipSnap)]);
disp(['totalSnapshots                  = ' num2str(totalSnapshots)]);
disp(['No of snapshots                 = ' num2str(nosSnapshot)]);

if PIconservative == 1 && (savePIvor == 1 || savePIt)
    disp('Conservative form of PI saved      ');
else
    disp('Non-Conservative form of PI saved  ');
end

disp(['No of snapshots                 = ' num2str(nosSnapshot)]);
disp(['DATA_DIR = ' DATA_DIR]);
disp('******************************************************');

%% DNS grid
N = 399;
NY = N + 1;
Lx = 6*pi;
dx = Lx/NX;
x = linspace(0,Lx-dx,NX);
kx = (2*pi/Lx)*[0:(NX/2) (-NX/2+1):-1];
Kx = ones(NY,1)*kx;
% kx(NX/2+1) = 0;
y = (-cos(pi*(0:N)/N))';
[X, Y] = meshgrid(x,y);

for countFilter = 1:length(filterTypeArr)
    filterType = convertStringsToChars(filterTypeArr(countFilter));

for count = 1:length(Ngrid)
    counter = 0;

    % Coarse grid
    NNX = Ngrid(count);
    ddx = Lx/NNX;
    xx = linspace(0,Lx-ddx,NNX);
    kkx = (2*pi/Lx)*[0:(NNX/2) (-NNX/2+1):-1];
    kkx(NNX/2+1) = 0;
    KKx = ones(NY,1)*kkx;
    [XX, YY] = meshgrid(xx,y);

    if saveVelocity == 1
        U = zeros(NY,NNX,totalSnapshots);
        V = zeros(NY,NNX,totalSnapshots);
    end
    
    if saveStreamfxnVorticity == 1
        Psi = zeros(NY,NNX,totalSnapshots);
        Vor = zeros(NY,NNX,totalSnapshots);
    end

    if saveTemperature == 1
        T = zeros(NY,NNX,totalSnapshots);
    end
    
    if saveResidualStress == 1
        S1 = zeros(NY,NNX,totalSnapshots);
        S2 = zeros(NY,NNX,totalSnapshots);
        S3 = zeros(NY,NNX,totalSnapshots);
    end
    
    if saveResidualStressVor == 1
        Svor1 = zeros(NY,NNX,totalSnapshots);
        Svor2 = zeros(NY,NNX,totalSnapshots);
    end

    if saveResidualStressT == 1
        J1 = zeros(NY,NNX,totalSnapshots);
        J2 = zeros(NY,NNX,totalSnapshots);
    end

    if savePIuv == 1
        PIu = zeros(NY,NNX,totalSnapshots);
        PIv = zeros(NY,NNX,totalSnapshots);
    end
    
    if savePIvor == 1
        PIvor = zeros(NY,NNX,totalSnapshots);
    end

    if savePIt == 1
        PIt = zeros(NY,NNX,totalSnapshots);
    end

    % Filtering
    
    if NNX == NX % LES grid is equal to DNS grid
        Delta = 2*Lx/Nfilter;
        N_LES = [NX NY];
    else
        Nfilter = Ngrid(count);
        Delta = 2*Lx/Nfilter;
        N_LES = [Ngrid(count) NY];
    end

    %%
    for countnosSnapshot = 1:length(nosSnapshot)
        
        % Code for Loading DNS data

        if saveVelocity == 1 || saveResidualStress == 1 || saveResidualStressVor == 1 ...
                || saveResidualStressT == 1 || savePIuv == 1 || savePIvor == 1
        
            load([DATA_DIR 'DNS' num2str(countnosSnapshot) '.mat'], 'slnStreamfunction');
            disp(['slnStreamfunction : ' DATA_DIR 'DNS' num2str(countnosSnapshot) '.mat loaded']);
        end
        
        if saveStreamfxnVorticity == 1 || saveResidualStressVor == 1 || savePIvor == 1
            load([DATA_DIR 'DNS' num2str(countnosSnapshot) '.mat'], 'slnVorticity');
            disp(['slnVorticity : ' DATA_DIR 'DNS' num2str(countnosSnapshot) '.mat loaded']);
        end

        if saveStreamfxnVorticity == 1 || saveResidualStressVor == 1 || savePIvor == 1
            load([DATA_DIR 'DNS' num2str(countnosSnapshot) '.mat'], 'slnTemperature');
            disp(['slnTemperature : ' DATA_DIR 'DNS' num2str(countnosSnapshot) '.mat loaded']);
        end
 
    for countSnap = 1:nosSnapshot(countnosSnapshot)   
        
        % Code to relate index of DNS data to filtered data
        
        if countnosSnapshot == 1
            counterTotal = countSnap;
        elseif countnosSnapshot > 1
            counterTotal = (countnosSnapshot-1)*nosSnapshot(1)+countSnap;
        end
        
        if mod(counterTotal,skipSnap) ~= 0
            counterTotal;
            continue;
        else
            counter = counter +1;
            [counterTotal counter];
        end

        % Derivatives
        order = 1; % Derivative order
        Dy = -cheb(N); % First derivative in y

        % Reshaping DNS data and initializing few DNS variables
        if saveVelocity == 1 || saveResidualStress == 1 || saveResidualStressVor == 1 ...
                || saveResidualStressT == 1 || savePIuv == 1 || savePIvor == 1

            Psi_DNS = reshape(slnStreamfunction(:,countSnap),NX,NY)';
            Psi_DNS_hat = fft(Psi_DNS,[],2);

            U_DNS_hat = zeros([NY,NX]);
            V_DNS_hat = zeros([NY,NX]);
            Ux_DNS_hat = zeros([NY,NX]);
            Vy_DNS_hat = zeros([NY,NX]);  

            U_DNS_hat(:,1) = (Dy)^order*Psi_DNS_hat(:,1);
            V_DNS_hat(:,1) = 0;
            Ux_DNS_hat(:,1) = 0;
            Vy_DNS_hat(:,1) = (Dy)^order*V_DNS_hat(:,1);
        end
        
        if saveStreamfxnVorticity == 1 || saveResidualStressVor == 1 || savePIvor == 1

            Vor_DNS = reshape(slnVorticity(:,countSnap),NX,NY)';
            Vor_DNS_hat = fft(Vor_DNS,[],2);

        end

        if savePIvor == 1 && PIconservative ~= 1
            % Non-conservative form of PI
            Vorx_DNS_hat = zeros([NY,NX]);
            Vory_DNS_hat = zeros([NY,NX]);

            Vorx_DNS_hat(:,1) = 0;
            Vory_DNS_hat(:,1) = (Dy)^order*Vor_DNS_hat(:,1);
        end

        if saveTemperature == 1 || savePIt == 1

            T_DNS = reshape(slnTemperature(:,countSnap),NX,NY)';
            T_DNS_hat = fft(T_DNS,[],2);
        end

        if savePIt == 1 && PIconservative ~= 1
            % Non-conservative form of PI
                Tx_DNS_hat = zeros([NY,NX]);
                Ty_DNS_hat = zeros([NY,NX]);

                Tx_DNS_hat(:,1) = 0;
                Ty_DNS_hat(:,1) = (Dy)^order*T_DNS_hat(:,1);
        end


        % Derivative in chebyshev grid
        for m = 2:NX/2
    
            if saveVelocity == 1 || saveResidualStress == 1 || saveResidualStressVor == 1 ...
                    || saveResidualStressT == 1 || savePIuv == 1 || savePIvor == 1

                U_DNS_hat(:,m) = (Dy)^order*Psi_DNS_hat(:,m);
                U_DNS_hat(:,NX-(m-2)) = conj(U_DNS_hat(:,m));
                
                V_DNS_hat(:,m) = -(1i*kx(m)).^order*Psi_DNS_hat(:,m);
                V_DNS_hat(:,NX-(m-2)) = conj(V_DNS_hat(:,m));

                Ux_DNS_hat(:,m) =(1i*kx(m)).^order*U_DNS_hat(:,m);
                Ux_DNS_hat(:,NX-(m-2)) = conj(Ux_DNS_hat(:,m));
                
                Vy_DNS_hat(:,m) = (Dy)^order*V_DNS_hat(:,m);
                Vy_DNS_hat(:,NX-(m-2)) = conj(Vy_DNS_hat(:,m));
       
            end
            
            if savePIvor == 1 && PIconservative ~= 1

                Vory_DNS_hat(:,m) = (Dy)^order*Vor_DNS_hat(:,m);
                Vory_DNS_hat(:,NX-(m-2)) = conj(Vory_DNS_hat(:,m));
                
                Vorx_DNS_hat(:,m) = (1i*kx(m)).^order*Vor_DNS_hat(:,m);
                Vorx_DNS_hat(:,NX-(m-2)) = conj(Vorx_DNS_hat(:,m));
    
            end
    
            if savePIt == 1 && PIconservative ~= 1
    
                Ty_DNS_hat(:,m) = (Dy)^order*T_DNS_hat(:,m);
                Ty_DNS_hat(:,NX-(m-2)) = conj(Ty_DNS_hat(:,m));
                
                Tx_DNS_hat(:,m) =(1i*kx(m)).^order*T_DNS_hat(:,m);
                Tx_DNS_hat(:,NX-(m-2)) = conj(Tx_DNS_hat(:,m));

            end
        end

        if saveVelocity == 1 || saveResidualStress == 1 || saveResidualStressVor == 1 ...
                || saveResidualStressT == 1 || savePIuv == 1 || savePIvor == 1

            U_DNS_hat(:,NX/2+1) = (Dy)^order*Psi_DNS_hat(:,NX/2+1);
            V_DNS_hat(:,NX/2+1) = -(1i*kx(NX/2+1)).^order*Psi_DNS_hat(:,NX/2+1);
            U_DNS = real(ifft(U_DNS_hat,[],2));
            V_DNS = real(ifft(V_DNS_hat,[],2));
    
            Ux_DNS_hat(:,NX/2+1) = (1i*kx(NX/2+1)).^order*U_DNS_hat(:,NX/2+1);
            Vy_DNS_hat(:,NX/2+1) = (Dy)^order*V_DNS_hat(:,NX/2+1);
            Ux_DNS = real(ifft(Ux_DNS_hat,[],2));
            Vy_DNS = real(ifft(Vy_DNS_hat,[],2));

            if NNX == NX
                U_filter = filter1D_RBC(U_DNS,filterType,'none', Delta, N_LES);
                V_filter = filter1D_RBC(V_DNS,filterType,'none', Delta, N_LES);
            else
                U_filter = filter1D_RBC(U_DNS,filterType,coarseGrainType, Delta, N_LES);
                V_filter = filter1D_RBC(V_DNS,filterType,coarseGrainType, Delta, N_LES);
            end
        end
        
        if saveStreamfxnVorticity == 1 || saveResidualStressVor == 1 || savePIvor == 1

            if NNX == NX
                Vor_filter = filter1D_RBC(Vor_DNS,filterType,'none', Delta, N_LES);
            else
                Vor_filter = filter1D_RBC(Vor_DNS,filterType,coarseGrainType, Delta, N_LES);
            end
        end

        if savePIvor == 1 && PIconservative ~= 1

            Vorx_DNS_hat(:,NX/2+1) = (1i*kx(NX/2+1)).^order*Vor_DNS_hat(:,NX/2+1);
            Vory_DNS_hat(:,NX/2+1) = (Dy)^order*Vor_DNS_hat(:,NX/2+1);
            Vorx_DNS = real(ifft(Vorx_DNS_hat,[],2));
            Vory_DNS = real(ifft(Vory_DNS_hat,[],2));

        end

        if saveTemperature == 1 || savePIt == 1

            if PIconservative ~= 1
                Tx_DNS_hat(:,NX/2+1) = (1i*kx(NX/2+1)).^order*T_DNS_hat(:,NX/2+1);
                Ty_DNS_hat(:,NX/2+1) = (Dy)^order*T_DNS_hat(:,NX/2+1);
                Tx_DNS = real(ifft(Tx_DNS_hat,[],2));
                Ty_DNS = real(ifft(Ty_DNS_hat,[],2));
            end

            if NNX == NX
                T_filter = filter1D_RBC(T_DNS,filterType,'none', Delta, N_LES);
            else
                T_filter = filter1D_RBC(T_DNS,filterType,coarseGrainType, Delta, N_LES);
            end
        end

        % Filtering and storing filtered Stream funxtion and Vorticity
        if saveStreamfxnVorticity == 1
            if NNX == NX
                Psi_filter = filter1D_RBC(Psi_DNS,filterType,'none', Delta, N_LES);
            else
                Psi_filter = filter1D_RBC(Psi_DNS,filterType,coarseGrainType, Delta, N_LES);
            end
            Psi(:,:,counter) = Psi_filter;
            Vor(:,:,counter) = Vor_filter;
        end

        % Filtering and storing filtered Velocity
        if saveVelocity == 1
            U(:,:,counter) = U_filter;
            V(:,:,counter) = V_filter;
        end

        if saveVelocity == 1
            U(:,:,counter) = U_filter;
            V(:,:,counter) = V_filter;
        end

        if saveTemperature == 1
            T(:,:,counter) = T_filter;
        end         

        if saveResidualStress == 1
            if NNX == NX
                UU_filter = filter1D_RBC(U_DNS.*U_DNS,filterType,'none', Delta, N_LES);
                VV_filter = filter1D_RBC(V_DNS.*V_DNS,filterType,'none', Delta, N_LES);
                UV_filter = filter1D_RBC(U_DNS.*V_DNS,filterType,'none', Delta, N_LES);
            else
                UU_filter = filter1D_RBC(U_DNS.*U_DNS,filterType,coarseGrainType, Delta, N_LES);
                VV_filter = filter1D_RBC(V_DNS.*V_DNS,filterType,coarseGrainType, Delta, N_LES);
                UV_filter = filter1D_RBC(U_DNS.*V_DNS,filterType,coarseGrainType, Delta, N_LES);
            end
            
            S1(:,:,counter) = UU_filter - U_filter.*U_filter;
            S2(:,:,counter) = UV_filter - U_filter.*V_filter;
            S3(:,:,counter) = VV_filter - V_filter.*V_filter;
        end  

        if saveResidualStressVor == 1
            if NNX == NX
                UVor_filter = filter1D_RBC(U_DNS.*Vor_DNS,filterType,'none', Delta, N_LES);
                VVor_filter = filter1D_RBC(V_DNS.*Vor_DNS,filterType,'none', Delta, N_LES);
            else
                UVor_filter = filter1D_RBC(U_DNS.*Vor_DNS,filterType,coarseGrainType, Delta, N_LES);
                VVor_filter = filter1D_RBC(V_DNS.*Vor_DNS,filterType,coarseGrainType, Delta, N_LES);
            end
            
            Svor1(:,:,counter) = UVor_filter - U_filter.*Vor_filter;
            Svor2(:,:,counter) = VVor_filter - V_filter.*Vor_filter;
        end

        if saveResidualStressT == 1
            if NNX == NX
                UT_filter = filter1D_RBC(U_DNS.*T_DNS,filterType,'none', Delta, N_LES);
                VT_filter = filter1D_RBC(V_DNS.*T_DNS,filterType,'none', Delta, N_LES);
            else
                UT_filter = filter1D_RBC(U_DNS.*T_DNS,filterType,coarseGrainType, Delta, N_LES);
                VT_filter = filter1D_RBC(V_DNS.*T_DNS,filterType,coarseGrainType, Delta, N_LES);
            end
            
            J1(:,:,counter) = UT_filter - U_filter.*T_filter;
            J2(:,:,counter) = VT_filter - V_filter.*T_filter;
        end      


%% For conservative PIvor and PIt

        if (savePIvor == 1 || savePIt == 1) && PIconservative == 1

            if savePIvor == 1

                UVorx_DNS_hat = zeros([NY,NX]);
                VVory_DNS_hat = zeros([NY,NX]);
                UVorx_filter_hat = zeros([NY,NNX]);
                VVory_filter_hat = zeros([NY,NNX]);
    
                UVor_DNS_hat  = fft(U_DNS.*Vor_DNS,[],2);
                VVor_DNS_hat  = fft(V_DNS.*Vor_DNS,[],2);
                UVor_filter_hat = fft(U_filter.*Vor_filter,[],2);
                VVor_filter_hat = fft(V_filter.*Vor_filter,[],2);

                UVorx_DNS_hat(:,1) = 0;
                VVory_DNS_hat(:,1) = (Dy)^order*VVor_DNS_hat(:,1);
            end

            if savePIt == 1

                UTx_DNS_hat = zeros([NY,NX]);
                VTy_DNS_hat = zeros([NY,NX]);
                UTx_filter_hat = zeros([NY,NNX]);
                VTy_filter_hat = zeros([NY,NNX]);
                
    
                UT_DNS_hat  = fft(U_DNS.*T_DNS,[],2);
                VT_DNS_hat  = fft(V_DNS.*T_DNS,[],2);
                UT_filter_hat = fft(U_filter.*T_filter,[],2);
                VT_filter_hat = fft(V_filter.*T_filter,[],2);

                UTx_DNS_hat(:,1) = 0;
                VTy_DNS_hat(:,1) = (Dy)^order*VT_DNS_hat(:,1);

            end
            
            for m = 2:NX/2

                if savePIvor == 1

                    VVory_DNS_hat(:,m) = (Dy)^order*VVor_DNS_hat(:,m);
                    VVory_DNS_hat(:,NX-(m-2)) = conj(VVory_DNS_hat(:,m));
                    
                    UVorx_DNS_hat(:,m) = (1i*kx(m)).^order*UVor_DNS_hat(:,m);
                    UVorx_DNS_hat(:,NX-(m-2)) = conj(UVorx_DNS_hat(:,m));
    
                end
                
                if savePIt == 1

                    VTy_DNS_hat(:,m) = (Dy)^order*VT_DNS_hat(:,m);
                    VTy_DNS_hat(:,NX-(m-2)) = conj(VTy_DNS_hat(:,m));
                    
                    UTx_DNS_hat(:,m) =(1i*kx(m)).^order*UT_DNS_hat(:,m);
                    UTx_DNS_hat(:,NX-(m-2)) = conj(UTx_DNS_hat(:,m));

                end
            end

            if savePIvor == 1

                UVorx_DNS_hat(:,NX/2+1) = (1i*kx(NX/2+1)).^order*UVor_DNS_hat(:,NX/2+1);
                VVory_DNS_hat(:,NX/2+1) = (Dy)^order*VVor_DNS_hat(:,NX/2+1);
                
                UVorx_filter_hat(:,1) = 0;
                VVory_filter_hat(:,1) = (Dy)^order*VVor_filter_hat(:,1);
            end

            if savePIt == 1

                UTx_DNS_hat(:,NX/2+1) = (1i*kx(NX/2+1)).^order*UT_DNS_hat(:,NX/2+1);
                VTy_DNS_hat(:,NX/2+1) = (Dy)^order*VT_DNS_hat(:,NX/2+1);
    
                UTx_filter_hat(:,1) = 0;
                VTy_filter_hat(:,1) = (Dy)^order*VT_filter_hat(:,1);
            end

            for m = 2:NNX/2

                if savePIvor == 1
                    VVory_filter_hat(:,m) = (Dy)^order*VVor_filter_hat(:,m);
                    VVory_filter_hat(:,NNX-(m-2)) = conj(VVory_filter_hat(:,m));
                
                    UVorx_filter_hat(:,m) = (1i*kkx(m)).^order*UVor_filter_hat(:,m);
                    UVorx_filter_hat(:,NNX-(m-2)) = conj(UVorx_filter_hat(:,m));

                end

                if savePIt == 1
                    
                    VTy_filter_hat(:,m) = (Dy)^order*VT_filter_hat(:,m);
                    VTy_filter_hat(:,NNX-(m-2)) = conj(VTy_filter_hat(:,m));
                    
                    UTx_filter_hat(:,m) =(1i*kkx(m)).^order*UT_filter_hat(:,m);
                    UTx_filter_hat(:,NNX-(m-2)) = conj(UTx_filter_hat(:,m));
                end
            end

            if savePIvor == 1
    
                UVorx_filter_hat(:,NNX/2+1) = (1i*kkx(NNX/2+1)).^order*UVor_filter_hat(:,NNX/2+1);
                VVory_filter_hat(:,NNX/2+1) = (Dy)^order*VVor_filter_hat(:,NNX/2+1);

                UVorx_DNS  = real(ifft(UVorx_DNS_hat,[],2));
                VVory_DNS  = real(ifft(VVory_DNS_hat,[],2));
    
                UVorx_filter = real(ifft(UVorx_filter_hat,[],2));
                VVory_filter = real(ifft(VVory_filter_hat,[],2));

            end

            if savePIt == 1

                UTx_filter_hat(:,NNX/2+1) = (1i*kkx(NNX/2+1)).^order*UT_filter_hat(:,NNX/2+1);
                VTy_filter_hat(:,NNX/2+1) = (Dy)^order*VT_filter_hat(:,NNX/2+1);
                
                UTx_DNS  = real(ifft(UTx_DNS_hat,[],2));
                VTy_DNS  = real(ifft(VTy_DNS_hat,[],2));
    
                UTx_filter = real(ifft(UTx_filter_hat,[],2));
                VTy_filter = real(ifft(VTy_filter_hat,[],2));

            end

        end

        if savePIvor == 1 

            if PIconservative == 1

                if NX == NNX
                    UVorx_VVory_filter_con = filter1D_RBC(UVorx_DNS + VVory_DNS,filterType,'none', Delta, N_LES);
                else
                    UVorx_VVory_filter_con = filter1D_RBC(UVorx_DNS + VVory_DNS,filterType,coarseGrainType, Delta, N_LES);
                end

                PIvor(:,:,counter) = UVorx_filter + VVory_filter - UVorx_VVory_filter_con;

            else
                if NX == NNX 
    
                    UVorx_filter_non = filter1D_RBC(U_DNS.*Vorx_DNS,filterType,'none', Delta, N_LES);
                    VVory_filter_non = filter1D_RBC(V_DNS.*Vory_DNS,filterType,'none', Delta, N_LES);
                    Vorx_filter = filter1D_RBC(Vorx_DNS,filterType,'none', Delta, N_LES);
                    Vory_filter = filter1D_RBC(Vory_DNS,filterType,'none', Delta, N_LES);
    
                else
    
                    UVorx_filter_non = filter1D_RBC(U_DNS.*Vorx_DNS,filterType,coarseGrainType, Delta, N_LES);
                    VVory_filter_non = filter1D_RBC(V_DNS.*Vory_DNS,filterType,coarseGrainType, Delta, N_LES);
                    Vorx_filter = filter1D_RBC(Vorx_DNS,filterType,coarseGrainType, Delta, N_LES);
                    Vory_filter = filter1D_RBC(Vory_DNS,filterType,coarseGrainType, Delta, N_LES);

                end

                PIvor(:,:,counter) = U_filter.*Vorx_filter + V_filter.*Vory_filter - UVorx_filter_non - VVory_filter_non;

            end
        end

        if savePIt == 1

            if PIconservative == 1 

                if NX == NNX

                    UTx_VTy_filter_con = filter1D_RBC(UTx_DNS + VTy_DNS,filterType,'none', Delta, N_LES);

                else
                    UTx_VTy_filter_con = filter1D_RBC(UTx_DNS + VTy_DNS,filterType,coarseGrainType, Delta, N_LES);
                end

                PIt(:,:,counter) = UTx_filter + VTy_filter - UTx_VTy_filter_con;

            else

                if NX == NNX

                    UTx_filter_non = filter1D_RBC(U_DNS.*Tx_DNS,filterType,'none', Delta, N_LES);
                    VTy_filter_non = filter1D_RBC(V_DNS.*Ty_DNS,filterType,'none', Delta, N_LES);
                    Tx_filter = filter1D_RBC(Tx_DNS,filterType,'none', Delta, N_LES);
                    Ty_filter = filter1D_RBC(Ty_DNS,filterType,'none', Delta, N_LES);

                else

                    UTx_filter_non = filter1D_RBC(U_DNS.*Tx_DNS,filterType,coarseGrainType, Delta, N_LES);
                    VTy_filter_non = filter1D_RBC(V_DNS.*Ty_DNS,filterType,coarseGrainType, Delta, N_LES);
                    Tx_filter = filter1D_RBC(Tx_DNS,filterType,coarseGrainType, Delta, N_LES);
                    Ty_filter = filter1D_RBC(Ty_DNS,filterType,coarseGrainType, Delta, N_LES);
                end

                PIt(:,:,counter) = U_filter .* Tx_filter + V_filter .* Ty_filter - UTx_filter_non - VTy_filter_non;

            end
        end
    end
    end

    if NNX == NX
        SAVE_DIR = ['../data/RBC/' dataType '/' filterType '/' ICtype '/NX' num2str(NNX) ...
            '/Nfilter' num2str(Nfilter)];
    else
        SAVE_DIR = ['../data/RBC/' dataType '/' filterType '/' ICtype '/NX' num2str(NNX)];
    end
    mkdir(SAVE_DIR);

    if saveVelocity == 1
        save([SAVE_DIR  '/UV.mat'], 'U', 'V', 'dummy', '-v7.3');
        size(U)
        size(V)
        disp([SAVE_DIR '/UV.mat saved']);
    end
    if saveStreamfxnVorticity == 1
        save([SAVE_DIR '/PsiVor.mat'], 'Psi', 'Vor', 'dummy', '-v7.3');
        size(Psi)
        size(Vor)
        disp([SAVE_DIR '/PsiVor.mat saved']);
    end
    if saveTemperature == 1
        save([SAVE_DIR '/T.mat'], 'T', 'dummy', '-v7.3');
        size(T)
        disp([SAVE_DIR '/T.mat saved']);
    end
    if saveResidualStress == 1
        save([SAVE_DIR '/S.mat'], 'S1', 'S2', 'S3', 'dummy', '-v7.3');
        size(S1)
        size(S2)
        size(S3)
        disp([SAVE_DIR '/S.mat saved']);
    end
    if saveResidualStressVor == 1
        save([SAVE_DIR '/Svor.mat'], 'Svor1', 'Svor2', 'dummy', '-v7.3');
        size(Svor1)
        size(Svor2)
        disp([SAVE_DIR '/Svor.mat saved']);
    end
    if saveResidualStressT == 1
        save([SAVE_DIR '/J.mat'], 'J1', 'J2', 'dummy', '-v7.3');
        size(J1)
        size(J2)
        disp([SAVE_DIR '/J.mat saved']);
    end
    if savePIuv == 1
        save([SAVE_DIR '/PIuv.mat'], 'PIu', 'PIv', 'dummy', '-v7.3');
        size(PIu)
        size(PIv)
        disp([SAVE_DIR '/PIuv.mat saved']);
    end
    if savePIvor == 1
        save([SAVE_DIR '/PIvor.mat'], 'PIvor', 'dummy', '-v7.3');
        size(PIvor)
        disp([SAVE_DIR '/PIvor.mat saved']);
    end
    if savePIt == 1
        save([SAVE_DIR '/PIt.mat'], 'PIt', 'dummy', '-v7.3');
        size(PIvor)
        disp([SAVE_DIR '/PIt.mat saved']);
    end
end
end