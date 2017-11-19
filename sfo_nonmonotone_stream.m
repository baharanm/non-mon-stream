% F: submodular function
% V: ground set
% C: weights
% k: upper bound on the cardinality of the desired set = rank of p-matchoid
% P: P in P-system constraint
% eps: error for estimating OPT
% beta = 1;
function [val, sset, eval_num, eval_num_par] = sfo_nonmonotone_stream(F, V, p_system, knapsacks, k, P, d, eps, eps_s,beta)
instances = [];
max_marginal = 0;
eval_num_par = 0;
best_element = [];

eval_num = 0;
for e = V
    %disp(e)
    flag = true;
    for p = 1:length(p_system)
        if ~satisfy_constraint(p_system{p}, e)
            flag = false;
            break;
        end
    end
    if ~flag || ~satisfy_constraint(knapsacks, e)
        continue
    end
    improve = inc(F,[],e);
    
    
    if (improve > max_marginal)
        max_marginal = improve;
        best_element = e;
        instances = update_buckets(F, instances, max_marginal, k, P, d, eps, eps_s);
    end
    
    max_eval = 0;
    for i = 1:numel(instances)
        [eval, instances(i)] = streaming_greedy_density(e, p_system, knapsacks, instances(i), 1, beta);
        max_eval = max(max_eval, eval);
        eval_num = eval_num + eval;
    end
    eval_num_par = eval_num_par + max_eval;
    
end

[val, sset, ~, eval] = get_best_solution(instances);
eval_num = eval_num + eval;
eval_num_par = eval_num_par + eval;
if (max_marginal > val)
    sset = best_element;
    val = max_marginal;
end
end

%% ================== functions =======================
function [eval_num, instance] = streaming_greedy_density(e, p_system, knapsacks, instance, buck_num, beta)
eval_num = 0;
F = instance.chain{buck_num};
S = get(F,'current_set');
[C, C_val, eval] = exchange_candidate(F, S, e, p_system, knapsacks);
eval_num = eval_num + eval;
new_set = [setdiff(S, C), e];
eval_num = eval_num + 1;
improve = inc(F, S, e) - get(F,'current_val');
%%%improve = inc(F, setdiff(S, C), e) - get(F,'current_val');%TODO: I changed this

flag = true;
for p = 1:length(p_system)
    if ~satisfy_constraint(p_system{p}, new_set)
        flag = false;
        break;
    end
end
if flag && satisfy_constraint(knapsacks, new_set) ...
        && improve / total_weight(knapsacks, e) >= instance.rho ...
        && improve >= instance.alpha + (1 + beta) * C_val 
    %new_set
    instance.chain{buck_num} = init(F, new_set);
    if isempty(C)
        return
    end
else
    C = e;
end    
for c = C
    if buck_num < length(instance.chain)
        [eval, instance] = streaming_greedy_density(c, p_system, knapsacks, instance, buck_num + 1, beta);
        eval_num = eval_num + eval;
    end
end
end


function [C, C_val, eval_num] = exchange_candidate(F, S, e, p_system, knapsacks)
C = [];
C_val = 0;
eval_num = 0;
for p = 1:length(p_system)
%for cnst = p_system
    cnst = p_system{p};
    min_improv = Inf;
    %%min_improv = -Inf;
    if (~satisfy_constraint(cnst, [S, e]))
        for i = 1:length(S)
            candidate = [S(1:i-1), S(i+1:end), e];
            if (satisfy_constraint(cnst, candidate) && satisfy_constraint(knapsacks, candidate))
                eval_num = eval_num + 1;
                improv = F(S(1:i)) - F(S(1:i-1));
                %%%improv = inc(F, candidate, []) - inc(F, [S(1:i-1), S(i+1:end)], []); %TODO: I changed this
                if (improv < min_improv)
                %%%if (improv > min_improv)
                    c_l = S(i);
                    min_improv = improv;
                end
            end
        end
        if min_improv < Inf && ~ismember(c_l, C)
        %%%if min_improv > -Inf && ~ismember(c_l, C)
            C = [C, c_l];
            C_val = C_val + min_improv; 
        end
    end
end
end

%% ==================  =======================
function [best_val, sset, instances, eval_num] = get_best_solution(instances)
best_val = 0;
sset = [];
eval_num = 0;
for i = 1:numel(instances)
    [instances(i), set, val, eval] = get_instance_val(instances(i));
    eval_num = eval_num + eval;
    if (val > best_val)
        best_val = val;
        sset = set;
    end
end
end


function [instance, best_set, best_val, eval_num] = get_instance_val(instance)
best_val = 0;
best_set = [];
eval_num = 0;
for i = 1:length(instance.chain)
    F = instance.chain{i};
    if isempty(get(F, 'current_set'))
        continue
    end
    [F_u, eval] = unconstrained_max(F);
    eval_num = eval_num + eval;
    val = get(F, 'current_val');
    val_u = get(F_u, 'current_val');
    if (val_u > val)
        instance.chain{i} = F_u;
        val = val_u;
        disp('======== unconstrained is better')
    end
    
    if (val > best_val)
        best_val = val;
        best_set = get(instance.chain{i}, 'current_set');
    end
end
end


function [Fx, eval] = unconstrained_max(F)
X = [];
eval = 0;
Y = get(F,'current_set');
if isempty(Y)
    return
end
Fx = init(F, X);
for i = length(Y):-1:1
    eval = eval + 2;
    a = inc(Fx, X, Y(i)) - get(Fx, 'current_val');
    b = F([Y(1:i-1), Y(i+1:end)]) - get(F, 'current_val');
    if (a >= b)
        X = [X, Y(i)];
        init(Fx, X);
    else
        Y(i) = [];
        init(F, Y);
    end
end
end

function instances = update_buckets(F, instances, max_marginal, k, P, d, eps, eps_s)
gamma = (2*(4*P-1)*max_marginal)/(4*P*(8*P+2*d-1));
R = gamma * (1+eps).^(0:log(k)/log(1+eps));
lambda = (eps_s * max_marginal) / (4 * k);
%A = lambda * (1+eps).^(0:log(2*k)/log(1+eps));
A = 2.^(ceil(log(lambda)):floor(log(lambda*2*k)));
R_new = R;
A_new = A;
R_old = [];

if (~isempty(instances))
    instances(extractfield(instances(:,1), 'rho') < R(1), :) = [];
end
if (~isempty(instances))
    instances(:, extractfield(instances(1,:), 'alpha') < A(1)) = [];
end
if (~isempty(instances))
    R_new(R <= instances(end,end).rho) = [];
    A_new(A <= instances(end,end).alpha) = [];
    R_old = extractfield(instances(:,1), 'rho');
end

[a, rn] = meshgrid(A, R_new);
[an, ro] = meshgrid(A_new, R_old);
chain = {repmat({init(F,[])}, P, 4*P)};
%chain = {repmat({init(F,[])}, 4*P, 1)};
%incs = {repmat({zeros(k, 1)}, 4*P, 1)};
arn_instances = struct('rho', num2cell(rn), 'alpha', num2cell(a), 'chain', repmat(chain, size(a)));%, 'incs', repmat(incs, size(a)));
anro_instances = struct('rho', num2cell(ro), 'alpha', num2cell(an), 'chain', repmat(chain, size(an)));%, 'incs', repmat(incs, size(an)));
instances = [[instances, anro_instances]; arn_instances];
end

