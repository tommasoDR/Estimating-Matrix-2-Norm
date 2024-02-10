% implementation of steepest gradient descent algorithm 

function [result, rel_gaps, vect_ng, time, iter] = SGD_Norm(A, x, epsilon, max_eval)
    
    % measure time for experiment
    tic;

    % fixed settings variable
    pause_iter = false;
    verbose = false;

    % initialize support variable
    rel_gaps = zeros(max_eval, 1);
    vect_ng = zeros(max_eval, 1);

    % check validity of input parameters
    check(A, x, epsilon, max_eval);
    
    % define function
    Q = transpose(A) * A;
    func = @(xTQx, xTx) -xTQx / xTx;
    
    % Compute norm of A with MATLAB library to have relative gap
    real_norm = norm(A);

    % define the gradient of the function
    gradient = @(x, fx, Qx, xTx) 2 * ((-fx) * x - Qx) / xTx;

    % Initialize variable
    iter = 0;
    alpha = 0;
    exit_status = "max_iter";

    while iter < max_eval
        
        % compute value of function, gradient and norm of gradient in x
        [Qx, xTQx, xTx] = compute_terms(x, Q);
        [norm_estimate, grad_values, norm_gradient] = evaluate(x, Qx, xTQx, xTx);
        vect_ng(iter+1) = norm_gradient;

        % compute relative gap
        rel_gap = (real_norm - norm_estimate) / real_norm;
        rel_gaps(iter+1) = abs(rel_gap);

        % print stats of current iteration
        if verbose
            fprintf( ['Iter %d - Rel_gap: %d\t Gradient norm: %d\t' ...
                'Step size: %d\n'], iter, rel_gap, norm_gradient, alpha);
        end

        if norm_gradient < epsilon
            exit_status = "optimal";
            iter = iter + 1;
            break;
        end

        % compute step size with exact line search
        alpha = compute_step_size(x, xTQx, xTx, grad_values);

        % move x
        x = x - alpha * grad_values;
     
        if pause_iter 
            pause;
        end

        iter = iter + 1;
        
    end
    
    if verbose & exit_status == "max_iter" 
        fprintf("\nThe maximum number of iterations has been reached. " + ...
            "Solution might not be accurate.\n")
    elseif verbose & exit_status == "optimal"
        fprintf("\nThe stopping criterion has been reached.\n")
    elseif verbose
        fprintf("\nSomething gone wrong...\n")
    end
    
    % results
    if exit_status == "max_iter" 
        iter = max_eval - 1;
    end
    
    result = norm_estimate;
    time = toc;
    rel_gaps = rel_gaps(1:iter);
    vect_ng = vect_ng(1:iter);
    
    

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
        error('Input max_eval must be a positive scalar value.');
    end
end

end