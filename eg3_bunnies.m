% 2021-04-13: test conformal on the bunnies

clc; clear; clf;
addpath(genpath('C:\Users\seems\Dropbox\ResearchProjects\MatlabToolbox\'));
%%
mesh_dir = 'C:\Users\seems\Dropbox\ResearchProjects\Dataset\bunny_conformal/';
all_meshes = dir([mesh_dir, '*.off']);
%% 20, 50
s1_id = 1;
s2_id = 218;
for s2_id = [20,50,100,150,200,210,218]
    
    S1 = MESH.MESH_IO.read_shape([mesh_dir, all_meshes(s1_id).name]);
    S2 = MESH.MESH_IO.read_shape([mesh_dir, all_meshes(s2_id).name]);
    
    
    S1 = MESH.compute_LaplacianBasis(S1, 100);
    S2 = MESH.compute_LaplacianBasis(S2, 100);
    
    
    figure(1); clf;
    subplot(1,2,1); plot_func_on_mesh(S1)
    subplot(1,2,2); plot_func_on_mesh(S2);
    
    
    para.num_samples = 2e3;
    S1.samples = MESH.get_samples(S1, para.num_samples, 'Euclidean');
    S2.samples = MESH.get_samples(S2, para.num_samples, 'Euclidean');
    
    B1 = S1.evecs(:,1:50);  B2 = S2.evecs(:, 1:50);
    rng(3);
    T21_ini = randi(S1.nv, S2.nv, 1);
    
    
    T21 = T21_ini;
    C12 = B2\B1(T21,:);
    T21 = knnsearch(B1(S1.samples,:)*C12', B2(S2.samples,:));
    for k = 2:50
        B11 = S1.evecs(S1.samples,1:k); B22 = S2.evecs(S2.samples, 1:k);
        Ev11 = S1.evals(1:k);  Ev22 = S2.evals(1:k);
        Ev11 = Ev11/sum(Ev11); Ev22 = Ev22/sum(Ev22); % normalize the delta to enforce the isometry
        for iter = 1:5
            C12 = B22\B11(T21,:);
            T21 = knnsearch(B11*diag(Ev11)*C12', B22*diag(Ev22));
            
            C12 = B22\B11(T21,:);
            T21 = knnsearch(B11*C12', B22);
        end
    end
    T21 = knnsearch(B1*C12', B2);
    
    figure(4);
    subplot(1,2,1); visualize_map_on_source(S2, S1, T21);
    subplot(1,2,2); visualize_map_on_target(S2, S1, T21);
    
    
    %%
    Ev1 = S1.evals(1:50); Ev2 = S2.evals(1:50);
    C12_ini = fMAP.pMap2fMap(B2, B1, T21_ini);
    x2 = optimize_fmap_unconstrained(C12_ini,'lap', Ev1, Ev2);
    C12= reshape(x2,50,50);
    T21_cont = fMAP.fMap2pMap(B1, B2, C12);
    C12_cont = fMAP.pMap2fMap(B2, B1, T21_cont);
    %%
    figure(5);
    subplot(1,4,1); visualize_map_on_source(S2, S1, T21);
    subplot(1,4,2); visualize_map_on_target(S2, S1, T21_ini);
    subplot(1,4,3); visualize_map_on_target(S2, S1, T21_cont);
    subplot(1,4,4); visualize_map_on_target(S2, S1, T21);
    
    save_dir = ['results/test0413_bunny/maps/']; if ~isdir(save_dir), mkdir(save_dir); end
    
    dlmwrite([save_dir, num2str(s2_id),'_ini.map'], T21_ini);
    dlmwrite([save_dir, num2str(s2_id),'_cont.map'], T21_cont);
    dlmwrite([save_dir, num2str(s2_id),'_ours.map'], T21);
    
end