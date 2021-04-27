function S = normalize_mesh_area(S_,a)
if nargin < 2, a = 1; end
S = S_;
S.surface.TRIV = S_.surface.TRIV;
S.surface.VERT = S_.surface.VERT/sqrt(sum(MESH.calc_tri_areas(S_.surface)))*sqrt(a); % rescale such that the mesh_area = 0.5
S.surface.X = S.surface.VERT(:,1);
S.surface.Y = S.surface.VERT(:,2);
S.surface.Z = S.surface.VERT(:,3);

end