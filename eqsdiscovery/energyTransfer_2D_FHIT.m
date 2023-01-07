function [P] = energyTransfer_2D_FHIT(U,V,S11,S12,S22)
%% Energy transfer of 2D_FHIT using SGS stress

    Ux = derivative_2D_FHIT(U,[1,0],'U');
    Uy = derivative_2D_FHIT(U,[0,1],'U');
    Vx = derivative_2D_FHIT(V,[1,0],'V');
    
    P = -(S11-S22).*Ux - S12.*(Uy+Vx);

end