% dimension
m = 1000;
n = 50;

% test parameters
test_iter = 100;
epsilon = 1e-8;
max_eval = 5000;

% init of support variable
i = 1;
gaps_SGD = zeros(test_iter,1);
iters_SGD = zeros(test_iter,1);
times_SGD = zeros(test_iter,1);

gaps_FR = zeros(test_iter,1);
iters_FR = zeros(test_iter,1);
times_FR = zeros(test_iter,1);

gaps_PR = zeros(test_iter,1);
iters_PR = zeros(test_iter,1);
times_PR = zeros(test_iter,1);

gaps_HS = zeros(test_iter,1);
iters_HS = zeros(test_iter,1);
times_HS = zeros(test_iter,1);


while i <= test_iter
    % random matrix vector
    A = sprand(m,n,0.05,0.5); A(A~=0) = 10*nonzeros(A); A = full(A);

    %A = rand(m,n) * 10 - 5;
    x = rand(n,1) * 10 - 5;
   
    % test
    [~, gap_SGD, ng_SGD, time_SGD, iter_SGD] = SGD_Norm(A, x, epsilon, max_eval);
    [~, gap_FR, ng_FR, time_FR, iter_FR] = CG_Norm(A, x, epsilon, max_eval, 1);
    [~, gap_PR, ng_PR, time_PR, iter_PR] = CG_Norm(A, x, epsilon, max_eval, 2);
    [~, gap_HS, ng_HS, time_HS, iter_HS] = CG_Norm(A, x, epsilon, max_eval, 3);
    
    % crunching results
    if gap_SGD(iter_SGD) == 0
        gap_SGD(iter_SGD) = 10e-17;
    end
    gaps_SGD(i) = floor(log10(abs(gap_SGD(iter_SGD))));
    times_SGD(i) = time_SGD;
    iters_SGD(i) = iter_SGD;

    if gap_FR(iter_FR) == 0
        gap_FR(iter_FR) = 1e-17;
    end
    gaps_FR(i) = floor(log10(abs(gap_FR(iter_FR))));
    times_FR(i) = time_FR;
    iters_FR(i) = iter_FR;

    if gap_PR(iter_PR) == 0
        gap_PR(iter_PR) = 1e-17;
    end
    gaps_PR(i) = floor(log10(abs(gap_PR(iter_PR))));
    times_PR(i) = time_PR;
    iters_PR(i) = iter_PR;

    if gap_HS(iter_HS) == 0
        gap_HS(iter_HS) = 1e-17;
    end
    gaps_HS(i) = floor(log10(abs(gap_HS(iter_HS))));
    times_HS(i) = time_HS;
    iters_HS(i) = iter_HS;

    i = i + 1;
end

% SGD
[mean_iter_SGD, max_iters_SGD, mean_gap_SGD, gaps_variance_SGD, ...
    max_gap_SGD, mean_time_SGD, times_variance_SGD, max_time_SGD] = compute_stats(iters_SGD, gaps_SGD, times_SGD);

% FR
[mean_iter_FR, max_iters_FR, mean_gap_FR, gaps_variance_FR, ...
    max_gap_FR, mean_time_FR, times_variance_FR, max_time_FR] = compute_stats(iters_FR, gaps_FR, times_FR);

% PR
[mean_iter_PR, max_iters_PR, mean_gap_PR, gaps_variance_PR, ...
    max_gap_PR, mean_time_PR, times_variance_PR, max_time_PR] = compute_stats(iters_PR, gaps_PR, times_PR);

% HS
[mean_iter_HS, max_iters_HS, mean_gap_HS, gaps_variance_HS, ...
    max_gap_HS, mean_time_HS, times_variance_HS, max_time_HS] = compute_stats(iters_HS, gaps_HS, times_HS);



function [mean_iter, max_iter, mean_gap, gaps_variance, max_gap,...
                mean_time, time_variance, max_time] = compute_stats(iters, gaps, times)
    mean_iter = mean(iters);
    max_iter = max(iters);

    mean_gap = mean(gaps);
    gaps_variance = var(gaps);
    max_gap = max(gaps);

    mean_time = mean(times);
    time_variance = var(times);
    max_time = max(times);
end
