function newScore = inc(F,A,el)
A = sfo_unique_fast(A);
F = init(F,A);

if sum(A==el)>0
    newScore = get(F,'current_val');
    return;
end
    Iv = diag([zeros(length(F.Y_loc), 1); ones(length(F.Gs_loc), 1)]);
    J = det(F.L_window([F.Y_loc; F.Gs_loc([A, el])], [F.Y_loc; F.Gs_loc([A, el])])) / det(F.L_window +  Iv);
newScore = J;
