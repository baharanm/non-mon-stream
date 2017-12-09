function [F,H] = init(F,sset)
sset = sfo_unique_fast(sset);
if ~isequal(sset,get(F,'current_set'))
    H = det(F.L(sset,sset));% / F.N;
    H = log(H);
    F = set(F,'current_val',H,'current_set',sset);
else
    H = get(F,'current_val');
end
