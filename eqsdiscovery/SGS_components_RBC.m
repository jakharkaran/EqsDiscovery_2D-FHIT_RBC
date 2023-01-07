function [J1_sgs, J1_leonard, J1_cross, J1_reynolds,...
    J2_sgs, J2_leonard, J2_cross, J2_reynolds] = SGS_components_RBC( ...
    U_DNS,V_DNS, T_DNS,filterType,coarseGrainingType, Delta, N_LES)

%Components of Residual Stress Term for RBC
% U_DNS = Uf (filtered/resolved) + Ud (residual)
% Uf_c filtered coarse grained
% Udf_c residual filtered coarse grained

% Output SGS, Leonard, cross, and Reynolds heat flux

Uf = filter1D_RBC(U_DNS,filterType,'none', Delta, N_LES);
Vf = filter1D_RBC(V_DNS,filterType,'none', Delta, N_LES);
Tf = filter1D_RBC(T_DNS,filterType,'none', Delta, N_LES);

Ud = U_DNS - Uf;
Vd = V_DNS - Vf;
Td = T_DNS - Tf;

Udf_c = filter1D_RBC(Ud,filterType,coarseGrainingType, Delta, N_LES);
Vdf_c = filter1D_RBC(Vd,filterType,coarseGrainingType, Delta, N_LES);
Tdf_c = filter1D_RBC(Td,filterType,coarseGrainingType, Delta, N_LES);

UdTd_f_c = filter1D_RBC(Ud.*Td,filterType,coarseGrainingType, Delta, N_LES);
VdTd_f_c = filter1D_RBC(Vd.*Td,filterType,coarseGrainingType, Delta, N_LES);

Uf_c = filter1D_RBC(U_DNS,filterType,coarseGrainingType, Delta, N_LES);
Vf_c = filter1D_RBC(V_DNS,filterType,coarseGrainingType, Delta, N_LES);
Tf_c = filter1D_RBC(T_DNS,filterType,coarseGrainingType, Delta, N_LES);

Uf_f_c = filter1D_RBC(Uf,filterType,coarseGrainingType, Delta, N_LES);
Vf_f_c = filter1D_RBC(Vf,filterType,coarseGrainingType, Delta, N_LES);
Tf_f_c = filter1D_RBC(Tf,filterType,coarseGrainingType, Delta, N_LES);

UT_f_c = filter1D_RBC(U_DNS.*T_DNS,filterType,coarseGrainingType, Delta, N_LES);
VT_f_c = filter1D_RBC(V_DNS.*T_DNS,filterType,coarseGrainingType, Delta, N_LES);

UfTf_f_c = filter1D_RBC(Uf.*Tf,filterType,coarseGrainingType, Delta, N_LES);
VfTf_f_c = filter1D_RBC(Vf.*Tf,filterType,coarseGrainingType, Delta, N_LES);

UfTd_f_c = filter1D_RBC(Uf.*Td,filterType,coarseGrainingType, Delta, N_LES);
UdTf_f_c = filter1D_RBC(Ud.*Tf,filterType,coarseGrainingType, Delta, N_LES);
VfTd_f_c = filter1D_RBC(Vf.*Td,filterType,coarseGrainingType, Delta, N_LES);
VdTf_f_c = filter1D_RBC(Vd.*Tf,filterType,coarseGrainingType, Delta, N_LES);

% SGS heat flux vector
J1_sgs = UT_f_c - Uf_c.*Tf_c;
J2_sgs = VT_f_c - Vf_c.*Tf_c;

% Leonard heat flux vector
J1_leonard = UfTf_f_c - Uf_f_c.*Tf_f_c;
J2_leonard = VfTf_f_c - Vf_f_c.*Tf_f_c;

% Cross heat flux vector
J1_cross = UfTd_f_c + UdTf_f_c - Uf_f_c.*Tdf_c - Udf_c.*Tf_f_c;
J2_cross = VfTd_f_c + VdTf_f_c - Vf_f_c.*Tdf_c - Vdf_c.*Tf_f_c;

% Reynolds heat flux vector
J1_reynolds = UdTd_f_c - Udf_c.*Tdf_c;
J2_reynolds = VdTd_f_c - Vdf_c.*Tdf_c;

end