function C = sfo_constraint_matroid(r,G,V,S)
C.G = G;
C.r = r;
C.V = V;
C.S = S;

C = class(C,'sfo_constraint_matroid', sfo_constraint);
end