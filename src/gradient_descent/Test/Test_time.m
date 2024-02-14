% Script that executed a time test with the specified parameters

% Matrix A dimensions
m = [2000, 10000];
n = [10, 100, 500, 1000];

% Test parameters
test_iter = 1;
epsilon = 1e-4;
max_eval = 5000;

% Initialization of support data structures
r = 1;
c = 1;
mean_times_SGD = zeros(numel(m), numel(n));
std_times_SGD = zeros(numel(m), numel(n));

mean_times_CG = zeros(numel(m), numel(n));
std_times_CG = zeros(numel(m), numel(n));

% First dimension
for d1 = m
    % Second dimension
    for d2 = n
        i = 1;
        times_SGD = zeros(test_iter,1);
        times_CG = zeros(test_iter,1);
      
        while i <= test_iter
            fprintf("Test num. %d for M = %d N = %d\n", i, d1, d2);

            % Random matrix and vector
            A = rand(d1, d2) * 20 - 10;
            x = rand(d2, 1) * 6 - 3;
           
            % Test
            [~, ~, ~, time_SGD, ~] = SGD_Norm(A, x, epsilon, max_eval);
            [~, ~, ~, time_CG, ~] = CG_Norm(A, x, epsilon, max_eval, 1);
            
            % crunching results
            times_SGD(i) = time_SGD;                     
            times_CG(i) = time_CG;       
        
            i = i + 1;
        end
        
        mean_times_SGD(r, c) = mean(times_SGD);
        std_times_SGD(r, c) = std(times_SGD);

        mean_times_CG(r, c) = mean(times_CG);
        std_times_CG(r, c) = std(times_CG);

        c = c +1;
    end

    c = 1;
    r = r + 1;
end

saveplot(mean_times_SGD(1,:), std_times_SGD(1,:), "SGD Time comparison - M = 2000", "SGD_2000", "c")
saveplot(mean_times_SGD(2,:), std_times_SGD(2,:), "SGD Time comparison - M = 10000", "SGD_10000", "r")

saveplot(mean_times_CG(1,:), std_times_CG(1,:), "CG Time comparison - M = 2000", "CG_2000", "c" )
saveplot(mean_times_CG(2,:), std_times_CG(2,:), "CG Time comparison - M = 10000", "CG_10000", "r")


% Saves the plot with mean and standard deviation of execution time
function [] = saveplot(mean, variance, title_plt, filename, color)
    bar(mean, color);
    hold on;
    errorbar(1:numel(mean), mean, variance, 'k.', 'LineWidth', 1); 
    xlabel('N');
    xticklabels({'10','100','500','1000'})
    ylabel('Time (s)');
    set(gca, 'YScale', 'log');
    title(title_plt);
    grid on;
    ax = gca;
    ax.YAxis.MinorTick = 'off';
    ax.YMinorGrid = 'off';
    saveas(ax, filename+".png", "png")
    hold off;
end
        
