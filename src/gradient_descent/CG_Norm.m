function [min_val] = CG_Norm(A, x, epsilon, max_eval)

    % check validity of input parameters
    check(A, x, epsilon, max_eval);

    % define function
    Q = transpose(A) * A;
    func = @(xTQx, xTx) -xTQx / xTx;

    % Compute norm of A with MATLAB library to have relative gap
    real_norm = norm(A);

    % define the gradient of the function
    gradient = @(x, fx, Qx, xTx) 2 * ((-fx) * x - Qx) / xTx;

    % Initialize iterations variable
    iter = 0;
    alpha = 0;
  
    while iter <= max_eval

        % compute value of function, gradient and norm of gradient in x
        [Qx, xTQx, xTx] = compute_terms(x, Q);
        [norm_estimate, grad_values, norm_gradient] = evaluate(x, Qx, xTQx, xTx);

        % compute relative gap
        rel_gap = (real_norm - norm_estimate) / real_norm;

        % print stats of current iteration
        fprintf( ['Iter %d - Rel_gap: %d\t Gradient norm: %d\t' ...
            'Step size: %d\n'], iter, rel_gap, norm_gradient, alpha);

        if norm_gradient < epsilon
            break;
        end

        % compute new direction 
        if iter == 0
            d = grad_values;
        else
            beta = (norm_gradient^2) / (norm_grad_pred^2);
            d = grad_values + beta * d_pred;
        end

        % compute step size with exact line search
        alpha = compute_step_size(x, xTQx, xTx, d);

        % move point
        x = x - alpha * d;

        % update variable for next iteration
        d_pred = d;
        norm_grad_pred = norm_gradient;

        %pause;
        iter = iter + 1;
    end
    
    min_val = norm_estimate;

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
    
    alpha =  min(r(r > 0));   
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