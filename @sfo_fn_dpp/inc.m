% Implementation by Andreas Krause (krausea@gmail.com)
% 
% Example: See sfo_fn.m and the tutorial script for more information

function newScore = inc(F,A,el)
A = sfo_unique_fast(A);
F = init(F,A);

if sum(A==el)>0
    newScore = get(F,'current_val');
    return;
end
    J = det(F.L([A,el],[A,el]));% / F.N;

newScore = log(J);
