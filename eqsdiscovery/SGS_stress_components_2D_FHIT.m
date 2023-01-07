function [S11_SGS_stress, S12_SGS_stress, S22_SGS_stress, ...
    S11_leonard, S12_leonard, S22_leonard, ...
    S11_cross, S12_cross, S22_cross,...
    S11_reynolds, S12_reynolds, S22_reynolds] = SGS_stress_components_2D_FHIT( ...
    U_DNS,V_DNS, filterType,coarseGrainingType, Delta, N_LES)

%Components of Residual Stress Term for 2D Turbulence
% U_DNS = Uf (filtered/resolved) + Ud (residual)
% Uf_c filtered coarse grained
% Udf_c residual filtered coarse grained

% Output SGS, Leonard, cross, and Reynolds stress

Uf = filter2D_2D_FHIT(U_DNS,filterType,'none', Delta, N_LES);
Vf = filter2D_2D_FHIT(V_DNS,filterType,'none', Delta, N_LES);

Ud = U_DNS - Uf;
Vd = V_DNS - Vf;

Udf_c = filter2D_2D_FHIT(Ud,filterType,coarseGrainingType, Delta, N_LES);
Vdf_c = filter2D_2D_FHIT(Vd,filterType,coarseGrainingType, Delta, N_LES);

UdUd_f_c = filter2D_2D_FHIT(Ud.*Ud,filterType,coarseGrainingType, Delta, N_LES);
VdVd_f_c = filter2D_2D_FHIT(Vd.*Vd,filterType,coarseGrainingType, Delta, N_LES);
UdVd_f_c = filter2D_2D_FHIT(Ud.*Vd,filterType,coarseGrainingType, Delta, N_LES);

Uf_c = filter2D_2D_FHIT(U_DNS,filterType,coarseGrainingType, Delta, N_LES);
Vf_c = filter2D_2D_FHIT(V_DNS,filterType,coarseGrainingType, Delta, N_LES);

Uf_f_c = filter2D_2D_FHIT(Uf,filterType,coarseGrainingType, Delta, N_LES);
Vf_f_c = filter2D_2D_FHIT(Vf,filterType,coarseGrainingType, Delta, N_LES);

UU_f_c = filter2D_2D_FHIT(U_DNS.*U_DNS,filterType,coarseGrainingType, Delta, N_LES);
VV_f_c = filter2D_2D_FHIT(V_DNS.*V_DNS,filterType,coarseGrainingType, Delta, N_LES);
UV_f_c = filter2D_2D_FHIT(U_DNS.*V_DNS,filterType,coarseGrainingType, Delta, N_LES);

UfUf_f_c = filter2D_2D_FHIT(Uf.*Uf,filterType,coarseGrainingType, Delta, N_LES);
VfVf_f_c = filter2D_2D_FHIT(Vf.*Vf,filterType,coarseGrainingType, Delta, N_LES);
UfVf_f_c = filter2D_2D_FHIT(Uf.*Vf,filterType,coarseGrainingType, Delta, N_LES);

UfUd_f_c = filter2D_2D_FHIT(Uf.*Ud,filterType,coarseGrainingType, Delta, N_LES);
VfVd_f_c = filter2D_2D_FHIT(Vf.*Vd,filterType,coarseGrainingType, Delta, N_LES);
UfVd_f_c = filter2D_2D_FHIT(Uf.*Vd,filterType,coarseGrainingType, Delta, N_LES);
VfUd_f_c = filter2D_2D_FHIT(Vf.*Ud,filterType,coarseGrainingType, Delta, N_LES);

% SGS stress tensor
S11_SGS_stress = UU_f_c - Uf_c.*Uf_c;
S12_SGS_stress = UV_f_c - Uf_c.*Vf_c;
S22_SGS_stress = VV_f_c - Vf_c.*Vf_c;

% Leonard stress tensor
S11_leonard = UfUf_f_c - Uf_f_c.*Uf_f_c;
S12_leonard = UfVf_f_c - Uf_f_c.*Vf_f_c;
S22_leonard = VfVf_f_c - Vf_f_c.*Vf_f_c;

% cross stress tensor
S11_cross = 2* (UfUd_f_c - Uf_f_c.*Udf_c);
S12_cross = UfVd_f_c + VfUd_f_c - Uf_f_c.*Vdf_c - Vf_f_c.*Udf_c;
S22_cross = 2* (VfVd_f_c - Vf_f_c.*Vdf_c);

% Reynolds stress tensor
S11_reynolds = UdUd_f_c - Udf_c.*Udf_c;
S12_reynolds = UdVd_f_c - Udf_c.*Vdf_c;
S22_reynolds = VdVd_f_c - Vdf_c.*Vdf_c;

end