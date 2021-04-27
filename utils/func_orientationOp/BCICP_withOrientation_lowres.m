function [T21_out, T12_out] = BCICP_withOrientation_lowres(S1,S2, B1, B2,...
    T21_highres, T12_highres, type, lowres_cache)
k1 = size(B1,2);
k2 = size(B2,2);
% simplify the mesh
[S1_lowres, T1_low2high, T1_high2low] = simplify_mesh_with_NNmap(S1,1e3,lowres_cache);
[S2_lowres, T2_low2high, T2_high2low] = simplify_mesh_with_NNmap(S2,1e3,lowres_cache);
% transfer the initial pMaps
T12 = T2_high2low(T12_highres(T1_low2high));
T21 = T1_high2low(T21_highres(T2_low2high));

for iter = 1:3
    % optimize with orientation-preservation term
    [T21, T12] = bcicp_refine_with_OrientationOp(S1_lowres,S2_lowres,...
        S1_lowres.evecs(:,1:k1),S2_lowres.evecs(:,1:k2),...
        T21,T12,1,type);
    % regular BCICP to improve the map: on simplified mesh
    [T21, T12] = bcicp_refine(S1_lowres,S2_lowres,...
        S1_lowres.evecs(:,1:k1),S2_lowres.evecs(:,1:k2),...
        T21,T12,2);
end
[T21, T12] = bcicp_refine(S1_lowres,S2_lowres,...
    S1_lowres.evecs(:,1:k1),S2_lowres.evecs(:,1:k2),...
    T21,T12,2);

% transfer to the original mesh
T12_new_highres = T2_low2high(T12(T1_high2low));
T21_new_highres = T1_low2high(T21(T2_high2low));
% BCICP refinement to improve the coverage
[T21_out, T12_out] = bcicp_refine(S1,S2,...
    B1,B2,T21_new_highres,T12_new_highres,2);
end