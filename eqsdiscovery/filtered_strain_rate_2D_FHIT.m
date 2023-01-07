function [SR11, SR12, SR22] = filtered_strain_rate_2D_FHIT(U,V)
%% Filtered rate of 2D_FHIT using SGS stress

% Input: filtered velocities

    Ux = derivative_2D_FHIT(U,[1,0],'U');
    Uy = derivative_2D_FHIT(U,[0,1],'U');
    Vx = derivative_2D_FHIT(V,[1,0],'V');

    SR11 = Ux;
    SR12 = 0.5*(Uy + Vx);
    SR22 = -Ux;

end