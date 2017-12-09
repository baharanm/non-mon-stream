function sat = satisfy_constraint(C, A)
sat = (sum(ismember(A, C.G)) == length(A));
end