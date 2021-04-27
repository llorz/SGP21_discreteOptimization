function [ SymOp, HOp ] = OrientationOp_local( surface, H)
%F Summary of this function goes here
%   Detailed explanation goes here

% compute its gradient
G = face_grads(surface, H);

% normalize it so that it has unit norm
vn  = sqrt(sum(G'.^2))';
vn = vn + eps;
vn = repmat(vn,1,3);
G = G./vn;

% rotate it by pi/2 using cross product with the normal
JGsrc = -J(surface,G);

% create 1st order differential operators associated with the vector fields
SymOp = surface.A*vf2op(surface, JGsrc);
HOp = surface.A*vf2op(surface, G);

end

