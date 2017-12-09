function sat = satisfy_constraint(C, A)
sat = sum(ismember([C.V(A)',C.S], C.G)) <= C.r;
end