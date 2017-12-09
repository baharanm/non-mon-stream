function C = sfo_constraint_knapsack(W,capacity)
C.W = W;
C.capacity = capacity;

C = class(C,'sfo_constraint_knapsack', sfo_constraint);
end