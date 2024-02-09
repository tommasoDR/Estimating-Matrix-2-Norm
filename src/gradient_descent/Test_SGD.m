% dimension
m = 1000;
n = 5;

% test parameters
test_iter = 1000;

% init of support variable
i = 1;
rel_gaps = zeros(test_iter,1);
iters = zeros(test_iter,1);

while i <= test_iter
    % random matrix and vector
    A = sprand(m,n,0.01,0.9);
    A(A~=0) = 10*nonzeros(A) - 5;
    A = full(A);
    %A = rand(m,n) * 10 - 5;
    x = rand(n,1) * 10 - 5;
    
    % norm gradient tolerance
    epsilon = 1e-6;
    
    % max evaluation
    max_eval = 5000;
    
    % test
    [min_value, rel_gap, iter] = SGD_Norm(A, x, epsilon, max_eval);
    rel_gaps(i) = floor(log10(abs(rel_gap)));

    iters(i) = iter;
    i = i + 1;
end

mean_iter = mean(iters)
max_iters = max(iters)
max_gap = max(rel_gaps)