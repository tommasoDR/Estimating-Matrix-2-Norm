function [min_val] = CG_Norm(A, x, epsilon, max_eval)

    % check validity of input parameters
    check(A, x, epsilon, max_eval);

    % define function
    Q = transpose(A) * A;
    func = @(x) -(x' * Q * x) / (x' * x);
    real_norm = norm(A);

    % define the gradient of the function
    gradient = @(x,fx) 2 * ((-fx) * x - Q * x) / (x' * x);
    
    % compute value of function, gradient and norm of gradient in x0
    [func_value, grad_values, norm_gradient] = evaluate(x);

    % Initialize stats variable
    func_eval = 1;
    iter = 1;
   
    init_flag = true;

    while norm_gradient >= epsilon && func_eval <= max_eval

        if init_flag
            d = grad_values;
            init_flag = false;
        else
            beta = (norm_gradient^2) / (norm_grad_pred^2);
            d = grad_values + beta * d_pred;
        end

        % compute stepsize 
        alpha = compute_step_size(x, d);

        % move point
        x = x - alpha * d;

        d_pred = d;
        norm_grad_pred = norm_gradient;

        rel_gap = (real_norm - sqrt(-func_value)) / real_norm;

        fprintf( 'Iter %d - Rel_gap: %d\t Func value: %d\t Gradient norm: %d\t Step size: %d\n', iter, rel_gap, func_value, norm_gradient, alpha);

        [func_value, grad_values, norm_gradient] = evaluate(x);

        %pause;
        iter = iter + 1;
        func_eval = func_eval + 1;
    end
    
    min_val = sqrt(-func_value);
    
% Compute step size
function [alpha] = compute_step_size(x,g)
    a = (g' * Q * g) * (g' * x) - (g' * Q * x) * (g' * g);
    b = (x' * Q * x) * (g' * g) - (g' * Q * g) * (x' * x);
    c = (x' * Q * g) * (x' * x) - (x' * Q * x) * (x' * g);
    
    r = roots([a b c]);

    alpha =  min(r(r > 0));
end

% Compute value of function, gradient and norm of gradient in x
function [func_value, grad_values, norm_gradient] = evaluate(x)
    func_value = func(x);
    grad_values = gradient(x, func_value);
    norm_gradient = norm(grad_values);
end


% Check input parameters
function [] = check(A,x0,epsilon,max_eval)

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