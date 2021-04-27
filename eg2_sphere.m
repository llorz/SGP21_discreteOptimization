clc; clear; clf
addpath(genpath('./utils/'))
addpath('func_main/')
mesh_dir = 'data/';

S1 = MESH.MESH_IO.read_shape([mesh_dir, 'sphere1.off']);
S1 = MESH.compute_LaplacianBasis(S1, 100);
S2 = S1;
S2.name = 'sphere2';
%% Input descriptors: we expect f1 to correspond to f2
vid1 = 20;
vid2 = 100;
% vid2 = 1260;
% vid2 = 939;

f1  = zeros(S1.nv,1); f1(vid1) = 1;
f2  = zeros(S2.nv,1); f2(vid2) = 1;

f1 = S1.evecs*pinv(S1.evecs)*f1;
f2 = S2.evecs*pinv(S2.evecs)*f2;

% visualize the two input descriptors
figure(1);
subplot(1,2,1); plot_func_on_mesh(S1, f1);
subplot(1,2,2); plot_func_on_mesh(S2, f2);
%% random initialization
rng(2);
T21_ini = randi(S1.nv, S2.nv,1);
%% Applying our discrete solver to minimize : |CC^T - I| + alpha*| Cf1 - f2|
all_T21 = {};
alpha = 1e3;
T21 = T21_ini;
all_T21{end+1} = T21;
for k = 2:50
    B1 = S1.evecs(:,1:k); B2 = S2.evecs(:, 1:k);
    f11 = pinv(B1)*f1; f22 = pinv(B2)*f2;
    for iter = 1:2
        C12 = B2\B1(T21,:);
        % the main update step
        T21 = knnsearch([B1*C12', alpha*B1*f11],...
            [B2, alpha*B2*f22]);
    end
    all_T21{end+1} = T21;
end

%% visualize the final map
figure(2);
subplot(1,2,1);
visualize_map_on_source(S2, S1, T21); title('Source')
subplot(1,2,2);
visualize_map_on_target(S2, S1, T21); title('Target')

%% visualize the mapped function 
f11 = pinv(B1)*f1;
g2 = B2*C12*f11;
figure(3);
subplot(1,3,1); plot_func_on_mesh(S1, f1); title('$f_1$','Interpreter','latex')
subplot(1,3,2); plot_func_on_mesh(S2, f2); title('$f_2$','Interpreter','latex')
subplot(1,3,3); plot_func_on_mesh(S2, g2); title('$C_{12}f_1$','Interpreter','latex')
%% visualize how f1 get transferred over iterations
all_C12 = cellfun(@(T21) fMAP.pMap2fMap(S2.evecs, S1.evecs, T21), all_T21,'uni', 0);
f11 = pinv(S1.evecs)*f1;
all_g = cellfun(@(C12) S2.evecs*C12*f11, all_C12, 'uni',0);
figure(4); clf;
check_id = [1:5,10,25,length(all_T21)];
for i = 1:length(check_id)
    subplot(2,4,i); plot_func_on_mesh(S2, all_g{check_id(i)}); title(['$C^k f_1 (k= $', num2str(check_id(i)),')'],'Interpreter','latex')
end
colormap(parula)
