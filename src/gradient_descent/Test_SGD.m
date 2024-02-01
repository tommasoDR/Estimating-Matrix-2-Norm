% dimension
m = 100;
n = 20;

% random matrix and vector
A = rand(m,n) * 4 - 2;
x = rand(n,1) * 10 - 5;

% norm gradient tolerance
epsilon = 1e-6;

% max evaluation
max_eval = 1000;

SGD_Norm(A, x, epsilon, max_eval)