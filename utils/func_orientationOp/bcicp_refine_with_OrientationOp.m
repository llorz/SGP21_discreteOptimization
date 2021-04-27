% 2018-05-10
function [T21, T12] = bcicp_refine_with_OrientationOp(S1,S2,B1,B2,T21_ini, T12_ini,num_iter,type)
compute_coverage = @(T) length(unique(T))/length(T);
T12 = T12_ini; T21 = T21_ini;

[S1.normals_vtx, S1.normals_face] = MESH.compute_vtx_and_face_normals(S1);
[S2.normals_vtx, S2.normals_face] = MESH.compute_vtx_and_face_normals(S2);

F11 = OrientationOp(S1,B1);
F22 = OrientationOp(S2,B2);

if compute_coverage(T12 < 0.5) && compute_coverage(T21 < 0.5)
    [T21, T12] = refine_pMap(T21,T12,S1,S2,4);
end
for iter = 1:num_iter
    C12 = B2\B1(T21,:);
    C21 = B1\B2(T12,:);
    % 2018-04-22 add OrientationOp
    C21 = optimize_fMap_with_OrientationOp(C21,B2,B1,F22,F11,type);
    C12 = optimize_fMap_with_OrientationOp(C12,B1,B2,F11,F22,type);    
    C12 = C21'*mat_projection(C21*C12);
    C21 = C12'*mat_projection(C12*C21);
    
    Y1 = [B1*C12',B1];
    Y2 = [B2, B2*C21'];
    
    d2 = dist_xy(Y1,Y2); % cost matrix S1 -> S2
    [~,T12] = min(d2,[],2);
    [~,T21] = min(d2',[],2);
        
    if compute_coverage(T12 < 0.5) && compute_coverage(T21 < 0.5)
        [T21, T12] = refine_pMap(T21,T12,S1,S2,4);
    end
    
    C1 = B1\B1(T21(T12),:);
    C1 = mat_projection(C1);
    C2 = B2\B2(T12(T21),:);
    C2 = mat_projection(C2);
    
    Y2 = [B1(T21,:)*C1',B2];
    Y1 = [B1, B2(T12,:)*C2'];
    d2 = dist_xy(Y2,Y1); % cost matrix S2 -> S1
    [~,T21] = min(d2,[],2);
    [~,T12] = min(d2',[],2);

    if compute_coverage(T12 < 0.5) && compute_coverage(T21 < 0.5)
        [T21, T12] = refine_pMap(T21,T12,S1,S2,4);
    end
end
end
%%
% 2018-05-10
% solve C21* = argmin orientation_term where
% orientation term: |C21 F22 - F11 C21 |_F^2
function [C21] = optimize_fMap_with_OrientationOp(C21,B2,B1,F22,F11,type)
assert(size(C21,1) == size(B1,2), 'Size does not match!');
assert(size(C21,2) == size(B2,2), 'Size does not match!');
k1 = size(C21,1); k2 = size(C21,2);

switch lower(type)
    case 'direct' %|C21 F22 - F11 C21|_F^2
        funObj = @(C21) deal(0.5*norm(reshape(C21,[k1,k2])*F22 - F11*reshape(C21,[k1,k2]),'fro')^2,...
            reshape(reshape(C21,[k1,k2])*(F22*F22') - F11'*reshape(C21,[k1,k2])*F22 - ...
            F11*reshape(C21,[k1,k2])*F22' + F11'*F11*reshape(C21,[k1,k2]),...
            [],1));
    case 'symmetric' % |C21 F22 + F11 C21|_F^2
        funObj = @(C21) deal(0.5*norm(reshape(C21,[k1,k2])*F22 + F11*reshape(C21,[k1,k2]),'fro')^2,...
            reshape(reshape(C21,[k1,k2])*(F22*F22') + F11'*reshape(C21,[k1,k2])*F22 + ...
            F11*reshape(C21,[k1,k2])*F22' + F11'*F11*reshape(C21,[k1,k2]),...
            [],1));
    otherwise
        error(['Unsupported energy type: ', type])
end

problem.options = optimoptions(@fminunc,...
    'Display','none',...
    'Algorithm','quasi-newton',...
    'SpecifyObjectiveGradient',true,...
    'MaxIterations',1e3);
problem.x0 = reshape(C21,[],1);
problem.objective = @(x) funObj(x);
problem.solver = 'fminunc';
C21 = reshape(fminunc(problem),[k1,k2]);
end
