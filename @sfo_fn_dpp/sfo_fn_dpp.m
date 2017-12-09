function F = sfo_fn_dpp(L)
F.L = L;
F.N = det(F.L + eye(size(F.L)));

F = class(F,'sfo_fn_dpp',sfo_fn);
F = set(F,'current_set',-1);
