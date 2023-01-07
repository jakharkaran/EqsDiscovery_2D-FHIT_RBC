function [Z] = enstrophyTransfer_2D_FHIT(Vor, Svor1, Svor2)
%% Enstrophy transfer of 2D_FHIT using SGS vorticity stress
    
    Vorx = derivative_2D_FHIT(Vor,[1,0],'Vor');
    Vory = derivative_2D_FHIT(Vor,[0,1],'Vor');
    
    Z = -Svor1.*Vorx - Svor2.*Vory;

end