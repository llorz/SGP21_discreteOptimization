clc; clear; clf
addpath(genpath('./utils/'))
addpath('func_main/')

mesh_dir = './data/SMAL/';
all_meshes = dir([mesh_dir, '*.obj']);

% the 50 shape pairs we reported in the paper
test_shapes = [ones(48,1), (2:49)';41, 3;41, 10];
iPair = 2;
%%
all_err = [];
s1_id = test_shapes(iPair, 1);
s2_id = test_shapes(iPair, 2);

S1 = MESH.MESH_IO.read_shape([mesh_dir, all_meshes(s1_id).name]);
S2 = MESH.MESH_IO.read_shape([mesh_dir, all_meshes(s2_id).name]);

S1 = normalize_mesh_area(S1, 1);
S2 = normalize_mesh_area(S2, 1);

S1 = MESH.compute_LaplacianBasis(S1, 100);
S2 = MESH.compute_LaplacianBasis(S2, 100);

rng(2);
T21_ini = randi(S1.nv, S2.nv, 1);
B1 = S1.evecs(:,1:50); B2 = S2.evecs(:,1:50);
Ev1 = S1.evals(1:50); Ev2 = S2.evals(1:50);

func_proj = @(C12)  fMAP.pMap2fMap(B2, B1, fMAP.fMap2pMap(B1, B2, C12));

figure(1);
subplot(1,2,1); plot_func_on_mesh(S1);
subplot(1,2,2); plot_func_on_mesh(S2);
%% E1:  Laplacian commutativity
E1 = @(C12) norm(C12*diag(Ev1) - diag(Ev2)*C12, 'fro');

C12_ini = fMAP.pMap2fMap(B2, B1, T21_ini);

% continuous solver
a = tic;
func = @(x) E1(reshape(x,50,50));
x0 = C12_ini(:);
options = optimoptions('fminunc','display','off');
x = fminunc(func, x0, options);
C12_cont = reshape(x, 50,50);
T21_cont = fMAP.fMap2pMap(B1, B2, C12_cont);
C12_cont_proj = func_proj(C12_cont);
t_cont = toc(a);
%% discrete solver
para.num_samples = 2e3;
S1.samples = MESH.get_samples(S1, para.num_samples, 'Euclidean');
S2.samples = MESH.get_samples(S2, para.num_samples, 'Euclidean');

a2 = tic;
T21 = T21_ini;
C12 = B2\B1(T21,:);
T21 = knnsearch(B1(S1.samples,:)*C12', B2(S2.samples,:));
for k = 4:50
    B11 = S1.evecs(S1.samples,1:k); B22 = S2.evecs(S2.samples, 1:k);
    Ev11 = S1.evals(1:k);  Ev22 = S2.evals(1:k);
    Ev11 = Ev11/sum(Ev11); Ev22 = Ev22/sum(Ev22); % normalize the delta to enforce the isometry
    for iter = 1:5
        C12 = B22\B11(T21,:);
        T21 = knnsearch(B11*diag(Ev11), B22*diag(Ev22)*C12);
        C12 = B22\B11(T21,:);
        T21 = knnsearch(B11*C12', B22);
    end
    
end
T21 = knnsearch(B1*C12', B2);
C12 = B2\B1(T21,:);
t_dis = toc(a2);
%%
figure(2);
subplot(1,4,1); visualize_map_on_source(S2, S1, T21_cont); title('Source')
subplot(1,4,2); visualize_map_on_target(S2, S1, T21_ini); title('randIni')
subplot(1,4,3); visualize_map_on_target(S2, S1, T21_cont); title('Continuous Solver')
subplot(1,4,4); visualize_map_on_target(S2, S1, T21); title('Ours')

[E1(C12_ini), E1(C12_cont_proj), E1(C12)]




