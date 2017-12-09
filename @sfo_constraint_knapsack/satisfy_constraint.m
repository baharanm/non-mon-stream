function sat = satisfy_constraint(C, A)
sat = sum(sum(C.W(A),1) <= C.capacity) == length(A);
end