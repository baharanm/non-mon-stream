function [F,H] = init(F,sset)
sset = sfo_unique_fast(sset);
if ~isequal(sset,get(F,'current_set'))
    Iv = diag([zeros(length(F.Y_loc), 1); ones(length(F.Gs_loc), 1)]);
    H = det(F.L_window([F.Y_loc; F.Gs_loc(sset)], [F.Y_loc; F.Gs_loc(sset)])) / det(F.L_window +  Iv);
    F = set(F,'current_val',H,'current_set',sset);
else
    H = get(F,'current_val');
end
