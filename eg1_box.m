clc; clear; clf
addpath(genpath('./utils/'))
addpath('func_main/')
mesh_dir = 'data/';


S1 = MESH.MESH_IO.read_shape([mesh_dir, 'box1.obj']);

% set type = 1: S1 and S2 have the same surface area
% set type = 2: S1 and S2 are conformal
type =1;
%% construct S2
if type == 1 % S1 and S2 have the same surface area
    S2 = S1;
    S2.surface.X = S2.surface.X/1.8;
    S2.surface.Y = S2.surface.Y/1.4;
    S2.surface.VERT = [S2.surface.X, S2.surface.Y, S2.surface.Z];
    
    
    S1 = normalize_mesh_area(S1, 1);
    S2 = normalize_mesh_area(S2, 1);
elseif type == 2 % S1 and S2 are conformal to each other
    S2 = S1;
    S1 = normalize_mesh_area(S1, 1);
    S2 = normalize_mesh_area(S2, 2);
end

S2.name = 'box2';

S1 = MESH.compute_LaplacianBasis(S1, 100);
S2 = MESH.compute_LaplacianBasis(S2, 100);

figure(1);
subplot(1,2,1); plot_func_on_mesh(S1);
subplot(1,2,2); plot_func_on_mesh(S2);
%% set random initialization
rng(2)
T21_ini = randi(S1.nv, S2.nv, 1);
%% Applying our discrete solver to optimize the area-perserving energy
% min | C C^T - I|
T21 = T21_ini;
for k = 2:100
    B11 = S1.evecs(:,1:k); B22 = S2.evecs(:, 1:k);
    for iter = 1:10
        C12 = B22\B11(T21,:);
        T21 = knnsearch(B11*C12', B22);
    end
end
T21_ortho = T21;
figure(3);
subplot(1,2,1); visualize_map_on_source(S2, S1, T21);
subplot(1,2,2); visualize_map_on_target(S2, S1, T21);
%% Applying our discrete solver to optimize the modified conformal energy
% min | C_12 D_1 C_12^T - D_2 |
T21 = T21_ini;
for k = 2:100
    B11 = S1.evecs(:,1:k); B22 = S2.evecs(:, 1:k);
    Ev11 = S1.evals(1:k);  Ev22 = S2.evals(1:k);
    Ev11 = Ev11/sum(Ev11); Ev22 = Ev22/sum(Ev22); 
    for iter = 1:10
        C12 = B22\B11(T21,:);
        T21 = knnsearch(B11*diag(Ev11)*C12', B22*diag(Ev22));
        
    end
end

figure(4);
subplot(1,2,1); visualize_map_on_source(S2, S1, T21);
subplot(1,2,2); visualize_map_on_target(S2, S1, T21);
T21_lap = T21;
%%  Applying the continuous solver  
B11 = S1.evecs(:,1:50); B22 = S2.evecs(:,1:50);
Ev1 = S1.evals(1:50); Ev2 = S2.evals(1:50);
C12_ini = fMAP.pMap2fMap(B22, B11, T21_ini);

% optimize the area-perserving energy
x1 = optimize_fmap_unconstrained(C12_ini,'ortho');
C12= reshape(x1,50,50);
T21_ortho_cont = fMAP.fMap2pMap(B11, B22, C12);


% optimize the standard conformal energy
x2 = optimize_fmap_unconstrained(C12_ini,'lap', Ev1, Ev2);
C12= reshape(x2,50,50);
T21_lap_cont = fMAP.fMap2pMap(B11, B22, C12);

%%
figure(5);
subplot(2,3,1); visualize_map_on_source(S2, S1, T21_ini); title('Source')
subplot(2,3,2); visualize_map_on_target(S2, S1, T21_ortho); title('Ours (area-preserving)')
subplot(2,3,3); visualize_map_on_target(S2, S1, T21_lap); title('Ours (conformal)')
subplot(2,3,4); visualize_map_on_target(S2, S1, T21_ini); title('randIni')
subplot(2,3,5); visualize_map_on_target(S2, S1, T21_ortho_cont); title('Continuous (area-presv.)')
subplot(2,3,6); visualize_map_on_target(S2, S1, T21_lap_cont); title('Continuous (conformal)')

