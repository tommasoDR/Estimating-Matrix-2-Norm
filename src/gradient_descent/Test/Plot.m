function [] = Plot(vector1, vector2, vector3, plt_title, filename)

    vector1(vector1==0) = vector1(vector1==0) + 1e-16;
    vector2(vector2==0) = vector2(vector2==0) + 1e-16;
    vector3(vector3==0) = vector3(vector3==0) + 1e-16;

    % Total number of iterations
    num_iterations1 = numel(vector1);
    num_iterations2 = numel(vector2);
    num_iterations3 = numel(vector3);
     
    % Create a vector representing the iterations
    iterations1 = 1:num_iterations1;
    iterations2 = 1:num_iterations2;
    iterations3 = 1:num_iterations3;
     
    % Plotting  
    if ~isempty(vector3)
        semilogy(iterations1, vector1, "r", iterations2, vector2, "g", iterations3, vector3, "b"); 
    else
        semilogy(iterations1, vector1, "r", iterations2, vector2, "b");
    end
   
    xlabel('Iterations');
    ylabel('Gradient norm');
    title(plt_title);
    
    if isempty(vector3)
        legend("SGD", "CG-FR");
    else
        legend("FR", "PR", "HS");
    end

    grid on
    ax = gca;
    ax.YAxis.MinorTick = 'off';
    ax.YMinorGrid = 'off';

    saveas(ax, filename+".png", "png")
end