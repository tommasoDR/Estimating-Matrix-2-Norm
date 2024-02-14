% Implementation of steepest gradient descent algorithm to estimate the
% norm 2 of a matrix A
%
% Input:
% - A: the matrix for which the norm 2 has to be estimated
% - x: the starting point
% - epsilon: gradient norm value for stopping criteria
% - max_eval: maximum number of function evaluation (= num. of iteration)
%
% Output:
% - results: estimate of norm(A,2)
% - rel_gaps: relative gap for each iteration
% - vect_gn: gradient norm for each iteration
% - time: execution time of the algorithm
% - iter: the number of iteration executed until a stopping condition has
%           been reached (convergence or max_eval)

function [result, rel_gaps, vect_gn, time, iter] = SGD_Norm(A, x, epsilon, max_eval)
    
    % Measuring time for experiment
    tic;

    % Fixed settings variable
    pause_iter = false;     % if true the iteration are interactive 
    verbose = false;         % if true at each iteration stats are printed

    % Initialize support data structures for experiments
    rel_gaps = zeros(max_eval, 1);
    vect_gn = zeros(max_eval, 1);

    % Check validity of input parameters
    check(A, x, epsilon, max_eval);
    
    % Define function
    Q = transpose(A) * A;
    func = @(xTQx, xTx) -xTQx / xTx;

    % Define the gradient of the function
    gradient = @(x, fx, Qx, xTx) 2 * ((-fx) * x - Qx) / xTx;
    
    % Compute norm of A with MATLAB library to calculate relative gap
    real_norm = norm(A);

    % Initialize support variables
    iter = 0;
    alpha = 0;
    exit_status = "max_iter";

    while iter < max_eval
        
        % Compute value of function, gradient and norm of gradient in x
        [Qx, xTQx, xTx] = compute_terms(x, Q);
        [norm_estimate, grad_values, norm_gradient] = evaluate(x, Qx, xTQx, xTx);

        % Save gradient norm at iteration iter for experiments
        vect_gn(iter+1) = norm_gradient;

        % Compute relative gap
        rel_gap = abs((real_norm - norm_estimate) / real_norm);

        % Save relative gap at iteration iter for experiments
        rel_gaps(iter+1) = rel_gap;

        % Print stats of current iteration
        if verbose
            fprintf( ['Iter %d - Rel_gap: %d\t Gradient norm: %d\t' ...
                'Step size: %d\n'], iter, rel_gap, norm_gradient, alpha);
        end

        % Stopping criteria
        if norm_gradient < epsilon
            exit_status = "optimal";
            iter = iter + 1;
            break;
        end

        % Compute step size with exact line search
        alpha = compute_step_size(x, xTQx, xTx, grad_values);

        % Move x
        x = x - alpha * grad_values;
     
        if pause_iter 
            pause;
        end

        iter = iter + 1;
        
    end
    
    if verbose & exit_status == "max_iter" 
        fprintf("\nThe maximum number of iterations has been reached. " + ...
            "Solution might not be accurate.\n\n");
    elseif verbose & exit_status == "optimal"
        fprintf("\nThe stopping criterion has been reached.\n\n");
    elseif verbose
        fprintf("\nSomething gone wrong...\n\n");
    end
    
    % Results
    if exit_status == "max_iter" 
        iter = max_eval - 1;
    end
    
    result = norm_estimate;
    time = toc;
    rel_gaps = rel_gaps(1:iter);
    vect_gn = vect_gn(1:iter);
    

% ---- Inner functions ----

% Compute the terms used in the succesive calculations
function [Qx, xTQx, xTx] = compute_terms(x,Q)
    Qx = Q * x;
    xTQx = x' * Qx;
    xTx = x' * x;
end

% Compute value of function, gradient and norm of gradient in x
function [norm_estimate, grad_values, norm_gradient] = evaluate(x, Qx, xTQx, xTx)
    func_value = func(xTQx,xTx);
    norm_estimate = sqrt(-func_value);
    grad_values = gradient(x, func_value, Qx, xTx);
    norm_gradient = norm(grad_values);
end

% Compute step size
function [alpha] = compute_step_size(x, xTQx, xTx, g)
    gTQg = (g' * Q * g);
    xTQg = (x' * Q * g);

    gTg = (g' * g);
    xTg = (x' * g);

    a = gTQg * xTg - xTQg * gTg;
    b = xTQx * gTg - gTQg * xTx;
    c = xTQg * xTx - xTQx * xTg;
    
    r = roots([a b c]);
    
    alpha =  min(r(r >= 0));  

end

% Check input parameters
function [] = check(A, x0, epsilon, max_eval)
    % Check if A is a matrix
    if ~ismatrix(A)
        error('A must be a matrix.');
    end

    % Check if x is a vector
    if ~isvector(x0) || ~isreal(x0)
        error('Starting x must be a vector.');
    end

    % Check if epsilon is a positive real value
    if ~isscalar(epsilon) || ~isreal(epsilon) || epsilon<=0
        error('Input epsilon must be a real positive scalar value.');
    end

    % Check if max_eval is a scalar value
    if  ~isreal(max_eval) || ~mod(max_eval,1) == 0 || max_eval<=0
        error('Input max_eval must be a positive integer value.');
    end
end

end