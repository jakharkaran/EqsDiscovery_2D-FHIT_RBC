addpath('../eqsdiscovery/');

%% Generate filtered and SGS flow fields for Two-Dimensional Forced Homogeneous Istropic Turbulence (2D-FHIT) on a square grid

% ICtype = 'train1'; 
% train1 train2 train3 train4 train5 test1

% dataType = 'Re20kNX1024nx0ny4r0p1';
% Re20kNX1024nx4ny0r0p1 Re20kNX1024nx4ny0r0p01 Re20kNX1024nx4ny4r0p1
% Re20kNX1024nx25ny25r0p1 Re100kNX2048nx4ny0r0p1

% Ngrid = [8,16,32,64,128];                    % Nfilter = Ngrid for coarse grid
% Nfilter = 128;     % Nfilter = Ngrid for coarse grid (NNX != NX)                            

% Select Filter type
% filterType = 'gaussian'; % gaussian box sharpSpectral gaussian+box
% coarseGrainType = 'spectral'; % spectral physical 

% Select filtered quantities to be saved | Equate to 1 to save 
saveVelocity = 1;   % Save filtered velocity (both U and V)
saveStreamfxnVorticity = 1;     % Save stream function and vorticity
saveResidualStress = 1;     % Save filtered residual stress \bar(uv) - \bar(u)\bar(v)
saveResidualStressVor = 1;    % Save filtered residual stress of vorticity \bar(uVor) - \bar(u)\bar(Vor)
savePIuv = 1;    % Save PIuv (Divergence of tau)
savePIvor = 1;    % Save PIvor (Curl of Divergence of tau)

% Specific to DNS data
% skipSnap = 1;  % Save every skipSnap th filtered snapshot

if strcmp(dataType,'Re20kNX1024nx4ny0r0p1') || strcmp(dataType,'Re20kNX1024nx4ny0r0p01') ...
        || strcmp(dataType,'Re20kNX1024nx4ny4r0p1') || strcmp(dataType,'Re20kNX1024nx25ny25r0p1')

    NX = 1024;
    nosSnapshot = ones(1,5)*2000; % number of snashots in each file used to save DNS data
    % There are 5 DNS files with 2000 snapshots each, total of 10000 snapshots

        if strcmp(ICtype,'test1')
        nosSnapshot = ones(1)*2000;
        end

elseif strcmp(dataType,'Re100kNX2048nx4ny0r0p1') || strcmp(dataType,'Re100kNX2048nx25ny25r0p1')

    NX = 2048;
    nosSnapshot = ones(1,20)*500;

    if strcmp(ICtype,'test1')
    nosSnapshot = ones(1,4)*500;
    end
end
totalSnapshots = floor(sum(nosSnapshot)/skipSnap);  % total number of filtered snapshots saved

DATA_DIR = ['../data/2D_FHIT/' dataType '/DNS/' ICtype '/'];

dummy = 'dummy';

disp('******************************************************');
disp('---------------------2D-FHIT-----------------------');
disp(['dataType                        = ' dataType]);
disp(['Filter Type                     = ' filterType]);
disp(['Coarse Grain Type               = ' coarseGrainType]);
disp(['ICtype                          = ' ICtype ]);
disp(['N DNS                           = ' num2str(NX)]);
disp(['Ngrid                           = ' num2str(Ngrid)]);
disp(['Nfilter                         = ' num2str(Nfilter)]);
disp(['skipSnap                        = ' num2str(skipSnap)]);
disp(['totalSnapshots                  = ' num2str(totalSnapshots)]);
disp(['No of snapshots                 = ' num2str(nosSnapshot)]);
disp(['DATA_DIR = ' DATA_DIR]);
disp('******************************************************');

%% DNS grid
Lx = 2*pi;
dx = Lx/NX;
x = linspace(0,Lx-dx,NX);
kx = (2*pi/Lx)*[0:(NX/2) (-NX/2+1):-1];

[Ky,Kx] = meshgrid(kx,kx);
Ksq = Kx.^2 + Ky.^2;
invKsq = 1./Ksq;
invKsq(1,1) = 0;
Kabs = sqrt(Ksq);

for count = 1:length(Ngrid)
    
    counter = 0;

    % Coarse grid
    NNX = Ngrid(count);
    kkx = (2*pi/Lx)*[0:(NNX/2) (-NNX/2+1):-1];
    [KKy,KKx] = meshgrid(kkx,kkx);
    KKsq = KKx.^2 + KKy.^2;
    
    if saveVelocity == 1
        U = zeros(NNX,NNX,totalSnapshots);
        V = zeros(NNX,NNX,totalSnapshots);
    end
    
    if saveStreamfxnVorticity == 1
        Psi = zeros(NNX,NNX,totalSnapshots);
        Vor = zeros(NNX,NNX,totalSnapshots);
    end
    
    if saveResidualStress == 1
        S1 = zeros(NNX,NNX,totalSnapshots);
        S2 = zeros(NNX,NNX,totalSnapshots);
        S3 = zeros(NNX,NNX,totalSnapshots);
    end
    
    if saveResidualStressVor == 1
        Svor1 = zeros(NNX,NNX,totalSnapshots);
        Svor2 = zeros(NNX,NNX,totalSnapshots);
    end
    
    if savePIuv == 1
        PIu = zeros(NNX,NNX,totalSnapshots);
        PIv = zeros(NNX,NNX,totalSnapshots);
    end
    
    if savePIvor == 1
        PIvor = zeros(NNX,NNX,totalSnapshots);
    end

    % Filtering
    
    if NNX == NX % FDNS grid is equal to DNS grid
        Delta = 2*Lx/Nfilter;
        N_LES = [NX NX];
    else
        Nfilter = Ngrid(count);
        Delta = 2*Lx/Nfilter;
        N_LES = [Ngrid(count) Ngrid(count)];
    end
    
    %%
    for countnosSnapshot = 1:length(nosSnapshot)
        
        % Code for Loading DNS data
        
        load([DATA_DIR 'DNS' num2str(countnosSnapshot) '.mat'], 'slnPsiDNS');
        disp(['slnPsiDNS : ' DATA_DIR 'DNS' num2str(countnosSnapshot) '.mat loaded']);
        
        if saveStreamfxnVorticity == 1 || saveResidualStressVor == 1 || savePIvor == 1
            load([DATA_DIR 'DNS' num2str(countnosSnapshot) '.mat'], 'slnWorDNS');
            disp(['slnWorDNS : ' DATA_DIR 'DNS' num2str(countnosSnapshot) '.mat loaded']);
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
            [counterTotal counter]
        end

        if saveVelocity == 1 || saveResidualStress == 1 || saveResidualStressVor == 1 ...
                || savePIvor == 1
        
            U_DNS = real(ifft2(1i*Ky.*fft2(slnPsiDNS(:,:,countSnap))));
            V_DNS = -real(ifft2(1i*Kx.*fft2(slnPsiDNS(:,:,countSnap))));
            
            if NNX == NX               
                U_filter = filter2D_2D_FHIT(U_DNS,filterType,'none', Delta, N_LES);
                V_filter = filter2D_2D_FHIT(V_DNS,filterType,'none', Delta, N_LES);                
            else           
                U_filter = filter2D_2D_FHIT(U_DNS,filterType,coarseGrainType, Delta, N_LES);
                V_filter = filter2D_2D_FHIT(V_DNS,filterType,coarseGrainType, Delta, N_LES);          
            end
        end
        
        if saveStreamfxnVorticity == 1 || saveResidualStressVor == 1 || savePIvor == 1
            if NNX == NX
                Vor_filter = filter2D_2D_FHIT(slnWorDNS(:,:,countSnap),filterType,'none', Delta, N_LES);
            else
                Vor_filter = filter2D_2D_FHIT(slnWorDNS(:,:,countSnap),filterType,coarseGrainType, Delta, N_LES);
            end
        end

        % Filtering and storing filtered Stream funxtion and Vorticity
        if saveStreamfxnVorticity == 1
            if NNX == NX
                Psi_filter = filter2D_2D_FHIT(slnPsiDNS(:,:,countSnap),filterType,'none', Delta, N_LES);
            else
                Psi_filter = filter2D_2D_FHIT(slnPsiDNS(:,:,countSnap),filterType,coarseGrainType, Delta, N_LES);
            end
            Psi(:,:,counter) = Psi_filter;
            Vor(:,:,counter) = Vor_filter;
        end
        
        % Filtering and storing filtered Velocity
        if saveVelocity == 1
            U(:,:,counter) = U_filter;
            V(:,:,counter) = V_filter;
        end
        
        if saveResidualStress == 1
            if NNX == NX
                UU_filter = filter2D_2D_FHIT(U_DNS.*U_DNS,filterType,'none', Delta, N_LES);
                VV_filter = filter2D_2D_FHIT(V_DNS.*V_DNS,filterType,'none', Delta, N_LES);
                UV_filter = filter2D_2D_FHIT(U_DNS.*V_DNS,filterType,'none', Delta, N_LES);
            else
                UU_filter = filter2D_2D_FHIT(U_DNS.*U_DNS,filterType,coarseGrainType, Delta, N_LES);
                VV_filter = filter2D_2D_FHIT(V_DNS.*V_DNS,filterType,coarseGrainType, Delta, N_LES);
                UV_filter = filter2D_2D_FHIT(U_DNS.*V_DNS,filterType,coarseGrainType, Delta, N_LES);
            end
            
            S1(:,:,counter) = UU_filter - U_filter.*U_filter;
            S2(:,:,counter) = UV_filter - U_filter.*V_filter;
            S3(:,:,counter) = VV_filter - V_filter.*V_filter;
        end   
        
        if saveResidualStressVor == 1
            if NNX == NX
                UVor_filter = filter2D_2D_FHIT(U_DNS.*slnWorDNS(:,:,countSnap),filterType,'none', Delta, N_LES);
                VVor_filter = filter2D_2D_FHIT(V_DNS.*slnWorDNS(:,:,countSnap),filterType,'none', Delta, N_LES);
            else
                UVor_filter = filter2D_2D_FHIT(U_DNS.*slnWorDNS(:,:,countSnap),filterType,coarseGrainType, Delta, N_LES);
                VVor_filter = filter2D_2D_FHIT(V_DNS.*slnWorDNS(:,:,countSnap),filterType,coarseGrainType, Delta, N_LES);
            end
            
            Svor1(:,:,counter) = UVor_filter - U_filter.*Vor_filter;
            Svor2(:,:,counter) = VVor_filter - V_filter.*Vor_filter;
        end
        
        if savePIuv == 1    
            if NNX == NX
                UUx_filter = filter2D_2D_FHIT(real(ifft2(1i*Kx.*fft2(U_DNS.*U_DNS))),filterType,'none', Delta, N_LES);
                UVy_filter = filter2D_2D_FHIT(real(ifft2(1i*Ky.*fft2(U_DNS.*V_DNS))),filterType,'none', Delta, N_LES);
                UVx_filter = filter2D_2D_FHIT(real(ifft2(1i*Kx.*fft2(U_DNS.*V_DNS))),filterType,'none', Delta, N_LES);
                VVy_filter = filter2D_2D_FHIT(real(ifft2(1i*Ky.*fft2(V_DNS.*V_DNS))),filterType,'none', Delta, N_LES);
            else
                UUx_filter = filter2D_2D_FHIT(real(ifft2(1i*Kx.*fft2(U_DNS.*U_DNS))),filterType,coarseGrainType, Delta, N_LES);
                UVy_filter = filter2D_2D_FHIT(real(ifft2(1i*Ky.*fft2(U_DNS.*V_DNS))),filterType,coarseGrainType, Delta, N_LES);
                UVx_filter = filter2D_2D_FHIT(real(ifft2(1i*Kx.*fft2(U_DNS.*V_DNS))),filterType,coarseGrainType, Delta, N_LES);
                VVy_filter = filter2D_2D_FHIT(real(ifft2(1i*Ky.*fft2(V_DNS.*V_DNS))),filterType,coarseGrainType, Delta, N_LES);
            end
            
            U_filter_U_filterx = real(ifft2(1i*KKx.*fft2(U_filter.*U_filter)));
            U_filter_V_filtery = real(ifft2(1i*KKy.*fft2(U_filter.*V_filter)));
            U_filter_V_filterx = real(ifft2(1i*KKx.*fft2(U_filter.*V_filter)));
            V_filter_V_filtery = real(ifft2(1i*KKy.*fft2(V_filter.*V_filter))); 
        
            PIu(:,:,counter) = UUx_filter + UVy_filter - (U_filter_U_filterx + U_filter_V_filtery);
            PIv(:,:,counter) = UVx_filter + VVy_filter - (U_filter_V_filterx + V_filter_V_filtery);
            
        end
        
        if savePIvor == 1
            
            conu1 = 1i*Kx.*fft2(U_DNS.*slnWorDNS(:,:,countSnap));
            conv1 = 1i*Ky.*fft2(V_DNS.*slnWorDNS(:,:,countSnap));
            
            if NNX == NX
                convectionF = filter2D_2D_FHIT(conu1+conv1,filterType,'none', Delta, N_LES);
            else
                convectionF = filter2D_2D_FHIT(conu1+conv1,filterType,coarseGrainType, Delta, N_LES);
            end
            
            conuF2 = 1i*KKx.*fft2(U_filter.*Vor_filter);
            convF2 = 1i*KKy.*fft2(V_filter.*Vor_filter);
            convectionF2 = real(ifft2(conuF2 + convF2));

            PIvor(:,:,counter) = convectionF - convectionF2;
        end
    end
    end
    
    % Saving data
    
    if NNX == NX
        SAVE_DIR = ['../data/2D_FHIT/' dataType '/' filterType '/' ICtype '/NX' num2str(NNX) ...
            '/Nfilter' num2str(Nfilter)];
    else
        SAVE_DIR = ['../data/2D_FHIT/' dataType '/' filterType '/' ICtype '/NX' num2str(NNX)];
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
end
