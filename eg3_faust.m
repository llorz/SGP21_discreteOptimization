% compare to the multiplicative operators
clc; clear; clf
addpath(genpath('./utils/'))
addpath('func_main/')
mesh_dir = 'data/FAUST/';

%%
k1 = 20; k2 = 20; numTimes = 50; skipSize = 50;
mesh_options = {'IfComputeLB',true,'numEigs',100,... % compute k LB basis
    'IfComputeNormals',true,... % compute vtx normals for orientation term
    'IfComputeGeoDist',false};  % do not compute the geodesic distance matrix
para.beta = 0;
%%
s1_name = 'tr_reg_000.off';
s2_name = 'tr_reg_044.off';

S1 = MESH.MESH_IO.read_shape([mesh_dir, s1_name]);
S1 = normalize_mesh_area(S1,1);
S1 = MESH.preprocess(S1, mesh_options{:});


S2 = MESH.MESH_IO.read_shape([mesh_dir, s2_name]);
S2 = normalize_mesh_area(S2,1);
S2 = MESH.preprocess(S2, mesh_options{:});
%% compute WKS descriptors
B1 = S1.evecs(:,1:k1); Ev1 = S1.evals(1:k1);
B2 = S2.evecs(:,1:k2); Ev2 = S2.evals(1:k2);

fct1 = waveKernelSignature(S1.evecs(:,1:100), S1.evals(1:100), S1.A, numTimes);
fct2 = waveKernelSignature(S2.evecs(:,1:100), S2.evals(1:100), S2.A, numTimes);

fct1 = fct1(:,1:skipSize:end);
fct2 = fct2(:,1:skipSize:end);

figure(1);clf;
subplot(1,2,1); plot_func_on_mesh(S1, fct1); view([0,90])
subplot(1,2,2); plot_func_on_mesh(S2, fct2); view([0,90])

%% multiplicative term 

para.beta = 0;
[C12_direct, ~, E1] = compute_fMap_regular_with_orientationOp_local(S1,S2,B1,B2,Ev1,Ev2,fct1,fct2,'direct',para);
T21_multi = fMAP.fMap2pMap(B1, B2, C12_direct);

%% Ours
alpha1 =  1e-1; alpha2 = 1e-1;
T21 = knnsearch(fct1, fct2);
Ev2_sum =   cumsum(S2.evals);
for k = 2:20
    B11 = S1.evecs(:,1:k);  Ev11 = S1.evals(1:k);
    B22 = S2.evecs(:,1:k);  Ev22 = S2.evals(1:k);
    Ev11 = Ev11/sum(Ev11); Ev22 = Ev22/sum(Ev22); % normalize the delta to enforce the isometry
    F1 = pinv(B11)*fct1; F2 = pinv(B22)*fct2;
    for iter = 1:5
        C12 = B22\B11(T21,:);
        % descriptor term + Laplacian term + ortho
        T21 = knnsearch( [alpha1*B11*F1,  alpha2*B11*diag(Ev11)*C12', B11*C12'], ....
                                     [alpha1*B22*F2,  alpha2*B22*diag(Ev22),           B22]);
    end
    
    C12 = B22\B11(T21,:);
    T21 = knnsearch(B11*C12', B22);
end
T21_ours= T21;
%%
figure(2); clf;
subplot(1,3,1); visualize_map_on_source(S2, S1, T21_multi); title('Source'); view([0,90])
subplot(1,3,2); visualize_map_on_target(S2, S1, T21_multi); title('Multiplicative Op'); view([0,90])
subplot(1,3,3); visualize_map_on_target(S2, S1, T21_ours); title('Ours'); view([0,90])

return
%% Continuous solver/Ours - with orientation term
para.beta = 1;
[C12_direct, para2, E2] = compute_fMap_regular_with_orientationOp_local(S1,S2,B1,B2,Ev1,Ev2,fct1,fct2,'direct',para);
T21_orient = fMAP.fMap2pMap(B1, B2, C12_direct);


Op1_full = OrientationOp_local(S1,fct1);
Op2_full = OrientationOp_local(S2,fct2);
% add the orientation operator
T21 = knnsearch(fct1, fct2);
alpha3 = para2.alpha;
for k = 2:20
    B11 = S1.evecs(:,1:k);  Ev11 = S1.evals(1:k);
    B22 = S2.evecs(:,1:k);  Ev22 = S2.evals(1:k);
    Ev11 = Ev11/sum(Ev11); Ev22 = Ev22/sum(Ev22); % normalize the delta to enforce the isometry
    % compute the orientation operator    
    alpha3 = para2.alpha*k/20;
    
    Op1 = B11'*S1.A*(repmat(fct1(:,1), [1,k]).*B11);
    Op2 = B22'*S2.A*(repmat(fct2(:,1), [1,k]).*B22);
    
    F1 = pinv(B11)*fct1; F2 = pinv(B22)*fct2;
    for iter = 1:5
        C12 = B22\B11(T21,:);
        % the third term is the orientation term
        T21 = knnsearch([alpha1*B11*F1, alpha2*B11*diag(Ev11)*C12', alpha3*B11*Op1*C12', B11*C12'], ...
                                    [alpha1*B22*F2, alpha2*B22*diag(Ev22),          alpha3*B22*Op2,          B22]);
    end
    
    C12 = B22\B11(T21,:);
    T21 = knnsearch(B11*C12', B22);
    
end
T21_orient_ours = T21;
%%
figure(3); clf;
subplot(1,3,1); visualize_map_on_source(S2, S1, T21_orient); title('Source'); view([0,90])
subplot(1,3,2); visualize_map_on_target(S2, S1, T21_orient_ours); title('Multiplicative Op'); view([0,90])
subplot(1,3,3); visualize_map_on_target(S2, S1, T21_orient_ours); title('Ours'); view([0,90])
