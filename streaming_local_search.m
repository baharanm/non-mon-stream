% Baharan Mirzasoleiman, Stefanie Jegelka, Andreas Krause, 
% "Streaming Non-monotone Submodular Maximization: Personalized Video 
% Summarization on the Fly", AAAI'18.


%% ================== Streaming Local Search =======================
% F: submodular function (we use the wrapper from Matlab Toolbox  
%    for Submodular Function Optimization (SFO).
%    SFO can be downloaded from https://las.inf.ethz.ch/sfo/)
% V: ground set.
% ind_systems: a collection of independence systems.
% knapsacks: multiple knapsack constraints.
% k: upper bound on the cardinality of the desired set
% P: p from a p-matchoid constraint (we use Streaming-Greedy from 
%    Chekuri et al. (2015) as Indstream.
% d: number of knapsack constraints.
% eps: error for estimating OPT.
% eps_s, beta: parameters of Streaming-Greedy, from Chekuri et al. (2015).
function [val, sset, eval_num, eval_num_par] = streaming_local_search(F, V, ind_systems, knapsacks, k, P, d, eps, eps_s, beta)
instances = [];
max_marginal = 0;
eval_num_par = 0;
best_element = [];

eval_num = 0;
for e = V
    flag = true;
    for p = 1:length(ind_systems)
        if ~satisfy_constraint(ind_systems{p}, e)
            flag = false;
            break;
        end
    end
    if ~flag || ~satisfy_constraint(knapsacks, e)
        continue
    end
    improve = inc(F,[],e);
    
    % lazily instantiate the thresholds. 
    if (improve > max_marginal)
        max_marginal = improve;
        best_element = e;
        instances = update_buckets(F, instances, max_marginal, k, P, d, eps, eps_s);
    end
    
    max_eval = 0;
    % parfor can be used for parallel processing.
    for i = 1:numel(instances)
        [eval, instances(i)] = streaming_greedy_density(e, ind_systems, knapsacks, instances(i), 1, beta);
        max_eval = max(max_eval, eval);
        eval_num = eval_num + eval;
    end
    eval_num_par = eval_num_par + max_eval;
    
end

% find the best solution among all instances of IndStream.
[val, sset, ~, eval] = get_best_solution(instances);
eval_num = eval_num + eval;
eval_num_par = eval_num_par + eval;
if (max_marginal > val)
    sset = best_element;
    val = max_marginal;
end
end

%% ================== Streaming Greedy with Density Threshold (IndStreamDensity) =======================
function [eval_num, instance] = streaming_greedy_density(e, ind_systems, knapsacks, instance, buck_num, beta)
eval_num = 0;
F = instance.chain{buck_num};
S = get(F,'current_set');
[C, C_val, eval] = exchange_candidate(F, S, e, ind_systems, knapsacks);
eval_num = eval_num + eval;
new_set = [setdiff(S, C), e];
eval_num = eval_num + 1;
improve = inc(F, setdiff(S, C), e) - get(F,'current_val');

flag = true;
for p = 1:length(ind_systems)
    if ~satisfy_constraint(ind_systems{p}, new_set)
        flag = false;
        break;
    end
end
if flag && satisfy_constraint(knapsacks, new_set) ...
        && improve / total_weight(knapsacks, e) >= instance.rho ...
        && improve >= instance.alpha + (1 + beta) * C_val 
    instance.chain{buck_num} = init(F, new_set);
    if isempty(C)
        return
    end
else
    C = e;
end    
% pass the discarded elements to the following instances of
% IndStreamDensity
for c = C
    if buck_num < length(instance.chain)
        [eval, instance] = streaming_greedy_density(c, ind_systems, knapsacks, instance, buck_num + 1, beta);
        eval_num = eval_num + eval;
    end
end
end

% Exchange candidate of Chekuri et al. (2015)
function [C, C_val, eval_num] = exchange_candidate(F, S, e, ind_systems, knapsacks)
C = [];
C_val = 0;
eval_num = 0;
for p = 1:length(ind_systems)
    cnst = ind_systems{p};
    min_improv = Inf;
    if (~satisfy_constraint(cnst, [S, e]))
        for i = 1:length(S)
            candidate = [S(1:i-1), S(i+1:end), e];
            if (satisfy_constraint(cnst, candidate) && satisfy_constraint(knapsacks, candidate))
                eval_num = eval_num + 1;
                improv = F(S(1:i)) - F(S(1:i-1));
                if (improv < min_improv)
                    c_l = S(i);
                    min_improv = improv;
                end
            end
        end
        if min_improv < Inf && ~ismember(c_l, C)
            C = [C, c_l];
            C_val = C_val + min_improv; 
        end
    end
end
end

%% ================== Functions =======================
% find the best solution among instances of IndStream
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

% Running unconstrained maximization on the solutions of all instances of IndStream
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
    end
    
    if (val > best_val)
        best_val = val;
        best_set = get(instance.chain{i}, 'current_set');
    end
end
end

% Unconstrained Submodular Maximization of Buchbinder et al. (2015)
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
% Density threshold
gamma = (2*max_marginal)/(4*P+4*sqrt(p)+d/sqrt(p)+1);
R = gamma * (1+eps).^(0:log(k)/log(1+eps));

% Estimating alpha, from Chekuri et al. (2015):
% A = {2^i: i \in Z} \cap [eps.m/4k, eps.m/2]
lambda = (eps_s * max_marginal) / (4 * k);
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
chain = {repmat({init(F,[])}, P, 2*ceil(sqrt(P))+1)};
arn_instances = struct('rho', num2cell(rn), 'alpha', num2cell(a), 'chain', repmat(chain, size(a)));
anro_instances = struct('rho', num2cell(ro), 'alpha', num2cell(an), 'chain', repmat(chain, size(an)));
instances = [[instances, anro_instances]; arn_instances];
end

