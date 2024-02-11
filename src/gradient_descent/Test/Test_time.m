% dimension
m = [1000, 10000];
n = [10, 100, 500, 1000];

% test parameters
test_iter = 100;
epsilon = 1e-7;
max_eval = 5000;

r = 1;
c = 1;
mean_times_SGD = zeros(numel(m), numel(n));
std_times_SGD = zeros(numel(m), numel(n));

mean_times_CG = zeros(numel(m), numel(n));
std_times_CG = zeros(numel(m), numel(n));

for d1 = m
    for d2 = n
        i = 1;
        times_SGD = zeros(test_iter,1);
        times_FR = zeros(test_iter,1);
        %times_PR = zeros(test_iter,1);
        %times_HS = zeros(test_iter,1);
        while i <= test_iter
            % random matrix vector       
            A = rand(d1,d2) * 20 - 10;
            x = rand(d2,1) * 6 - 3;
           
            % test
            [~, ~, ~, time_SGD, ~] = SGD_Norm(A, x, epsilon, max_eval);
            [~, ~, ~, time_FR, ~] = CG_Norm(A, x, epsilon, max_eval, 1);
            %[~, ~, ~, time_PR, ~] = CG_Norm(A, x, epsilon, max_eval, 2);
            %[~, ~, ~, time_HS, ~] = CG_Norm(A, x, epsilon, max_eval, 3);
            
            % crunching results
            
            times_SGD(i) = time_SGD;          
           
            times_FR(i) = time_FR;
            
            %times_PR(i) = time_PR;

            %times_HS(i) = time_HS;
            
        
            i = i + 1
        end
        
        mean_times_SGD(r,c) = mean(times_SGD);
        std_times_SGD(r,c) = std(times_SGD);

        mean_times_CG(r,c) = mean(times_FR);
        std_times_CG(r,c) = std(times_FR);

        %mean_time_PR = mean(times_PR);
        %time_std_PR = std(times_PR);

        %mean_time_HS = mean(times_HS);
        %time_std_HS = std(times_HS);

        c = c +1;
    end

    c = 1;
    r = r + 1;
end

saveplot(mean_times_SGD(1,:), std_times_SGD(1,:), "SGD Time comparison - M = 1000", "SGD_1000", "c")
saveplot(mean_times_SGD(2,:), std_times_SGD(2,:), "SGD Time comparison - M = 10000", "SGD_10000", "r")

saveplot(mean_times_CG(1,:), std_times_CG(1,:), "CG Time comparison - M = 1000", "CG_1000", "c" )
saveplot(mean_times_CG(2,:), std_times_CG(2,:), "CG Time comparison - M = 10000", "CG_10000", "r")

function [] = saveplot(mean, variance, title_plt, filename, color)
    figure;
    bar(mean, color);
    hold on;
    errorbar(1:numel(mean), mean, variance, 'k.', 'LineWidth', 1); 
    xlabel('N');
    xticklabels({'10','100','500','1000'})
    ylabel('Time (s)');
    set(gca, 'YScale', 'log');
    title(title_plt);
    %legend('Mean', 'Standard deviation', 'Location', 'Top-left');
    grid on;
    ax = gca;
    ax.YAxis.MinorTick = 'off';
    ax.YMinorGrid = 'off';
    saveas(ax, filename+".png", "png")
    hold off;
end
        
