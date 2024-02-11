% dimension
m = 1000;
n = 1000;

% test parameters
test_iter = 1;
epsilon = 1e-9;
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
    %A = sprand(m,n,1,1e-5);
    %A = full(A);
    %A(A~=0) = 10 * nonzeros(A);

    A = rand(m,n) * 20 - 10;
    x = rand(n,1) * 6 - 3;
   
    % test
    [~, gap_SGD, ng_SGD, time_SGD, iter_SGD] = SGD_Norm(A, x, epsilon, max_eval);
    [~, gap_FR, ng_FR, time_FR, iter_FR] = CG_Norm(A, x, epsilon, max_eval, 1);
    [~, gap_PR, ng_PR, time_PR, iter_PR] = CG_Norm(A, x, epsilon, max_eval, 2);
    [~, gap_HS, ng_HS, time_HS, iter_HS] = CG_Norm(A, x, epsilon, max_eval, 3);
    
    % crunching results
    if gap_SGD(iter_SGD) == 0
        gap_SGD(iter_SGD) = 10e-16;
    end
    gaps_SGD(i) = floor(log10(abs(gap_SGD(iter_SGD))));
    times_SGD(i) = time_SGD;
    iters_SGD(i) = iter_SGD;

    if gap_FR(iter_FR) == 0
        gap_FR(iter_FR) = 1e-16;
    end
    gaps_FR(i) = floor(log10(abs(gap_FR(iter_FR))));
    times_FR(i) = time_FR;
    iters_FR(i) = iter_FR;

    if gap_PR(iter_PR) == 0
        gap_PR(iter_PR) = 1e-16;
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

    i = i + 1
end

% SGD
[mean_iter_SGD, iters_std_SGD, mean_gap_SGD, gaps_std_SGD, ...
    mean_time_SGD, times_std_SGD] = compute_stats(iters_SGD, gaps_SGD, times_SGD);

% FR
[mean_iter_FR, iters_std_FR, mean_gap_FR, gaps_std_FR, ...
    mean_time_FR, times_std_FR] = compute_stats(iters_FR, gaps_FR, times_FR);

% PR
[mean_iter_PR, iters_std_PR, mean_gap_PR, gaps_std_PR, ...
    mean_time_PR, times_std_PR] = compute_stats(iters_PR, gaps_PR, times_PR);

% HS
[mean_iter_HS, iters_std_HS, mean_gap_HS, gaps_std_HS, ...
    mean_time_HS, times_std_HS] = compute_stats(iters_HS, gaps_HS, times_HS);



function [mean_iter, iters_std, mean_gap, gaps_std,...
                mean_time, time_std] = compute_stats(iters, gaps, times)
    mean_iter = mean(iters);
    iters_std = std(iters);

    mean_gap = mean(gaps);
    gaps_std = std(gaps);

    mean_time = mean(times);
    time_std = std(times);
end
