% implementation of steepest gradient descent algorithm 

function [min_val] = SGD_Norm(A, x, epsilon, max_eval, alpha_start, ...
    armijo_m, wolfe_m, tau)
   
    % check validity of input parameters
    check(A, x, epsilon, max_eval, alpha_start, armijo_m, wolfe_m, tau);
    
    % define function
    Q = transpose(A) * A;
    %func = @(x) x'*Q*x;
    func = @(x) -(x' * Q * x) / (x' * x);

    % define the gradient of the function
    %gradient = @(x,fx) 2*Q*x;
    gradient = @(x,fx) 2 * (fx * x - Q * x) / (x' * x);
    
    % compute value of function, gradient and norm of gradient in x0
    [func_value, grad_values, norm_gradient] = evaluate(x);

    % Initialize stats variable
    func_eval = 1;
    iter = 1;

    while norm_gradient > epsilon && func_eval <= max_eval
        
        % computing efficiently derivative of tomography in 0 with gradient as direction
        phi_deriv_zero = - (norm_gradient ^ 2); 

        % compute stepsize with AW or backtracking
        if wolfe_m == 0
            alpha = backtracking_linesearch(func_value, phi_deriv_zero, armijo_m, alpha_start, tau);
        else
            alpha = AW_linesearch(func_value, phi_deriv_zero, armijo_m, wolfe_m, alpha_start, tau);
        end

        x = x - alpha * grad_values;

        [func_value, grad_values, norm_gradient] = evaluate(x);

        fprintf( 'Iter %d - Func eval: %d\t Func value: %d\t Gradient norm: %d\t Step size: %d\n', iter, func_eval, func_value, norm_gradient, alpha);
        
        %pause;
        func_eval = func_eval + 1;
        iter = iter + 1;
    end
    
    min_val = sqrt(-func_value);



% Backtraking line search
function [alpha] = backtracking_linesearch(func_value, phi_der_zero, armijo_m, alpha, tau)
    phi_alpha = func(x - alpha * grad_values);
    func_eval = func_eval + 1;

    while phi_alpha > (func_value + armijo_m * alpha * phi_der_zero)
        phi_alpha = func(x - alpha * grad_values);
        alpha = tau * alpha;
        func_eval = func_eval + 1;
    end
end

function [alpha] = AW_linesearch(func_value, phi_der_zero, armijo_m, wolfe_m, alpha, tau)
    while func_eval <= max_eval
        % evaluate tomography
        phi_alpha = func(x - alpha * grad_values);
        phi_der_alpha = - grad_values' * gradient(x - alpha * grad_values, phi_alpha);
        func_eval = func_eval + 1;

        % Armijo not satisfied
        if phi_alpha > (func_value + armijo_m * alpha * phi_der_zero)
            break;
        end

        % Wolfe satisfied
        if phi_der_alpha >= wolfe_m * phi_der_zero
            %alpha = alpha;
            return;
        end

        if phi_der_alpha >= 0
            break;
        end

        % Increment alpha
        alpha = alpha / tau;
   
    end
    
    min_alpha = 0;
    max_alpha = alpha;
    phi_der_min = phi_der_zero;
    phi_der_max = phi_der_alpha;

  
    while func_eval <= max_eval
        if phi_der_min < 0 && phi_der_max > 0
            % Quadratic interpolation
            fprintf("OKOK")
            alpha = (min_alpha * phi_der_max - max_alpha * phi_der_min) / (phi_der_max - phi_der_min);
            alpha = max([min_alpha + 0.05 * (max_alpha - min_alpha)  min([max_alpha - 0.05 * (max_alpha - min_alpha)  alpha])]);
        else
            % Binary search if quadratic interpolation is not applicable
            alpha = (min_alpha + max_alpha) / 2;
        end
        
        % evaluate tomography
        phi_alpha = func(x - alpha * grad_values);
        phi_der_alpha = - grad_values' * gradient(x - alpha * grad_values, phi_alpha);
        func_eval = func_eval + 1;

        if phi_alpha <= (func_value + armijo_m * alpha * phi_der_zero)
            if phi_der_alpha >= wolfe_m * phi_der_zero
                return;
            end

            min_alpha = alpha;
            phi_der_min = phi_der_alpha;
        else
            max_alpha = alpha;
            phi_der_max = phi_der_alpha;
        end

    end

end


% Compute value of function, gradient and norm of gradient in x
function [func_value, grad_values, norm_gradient] = evaluate(x)
    func_value = func(x);
    grad_values = gradient(x, func_value);
    norm_gradient = norm(grad_values);
end


% Check input parameters
function [] = check(A,x0,epsilon,max_eval,alpha_start, ...
    armijo_m,wolfe_m,tau)
    % Check if f is a function handle
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
    
    % Check if alpha_start is a positive real value
    if ~isscalar(alpha_start) || ~isreal(alpha_start) || alpha_start<=0
        error('Input alpha_start must be a real positive scalar value.');
    end

    % Check if armijo_m is a real value 
    if ~isscalar(armijo_m) || ~isreal(armijo_m) || armijo_m<=0 || armijo_m>=1
        error('Input armijo_m must be a real positive scalar value.');
    end

    % Check if wolfe_m is a real value 
    if ~isscalar(wolfe_m) || ~isreal(wolfe_m) || wolfe_m < 0 || wolfe_m>=1
        error('Input wolfe_m must be a real positive scalar value.');
    end
    
    %check values for armijo_m and wolfe_m
    if wolfe_m ~= 0 && armijo_m >= wolfe_m
        error('armijo_m and wolfe_m not valid')
    end
    
    % Check if tau is a real value 
    if ~isscalar(tau) || ~isreal(tau) || tau<=0 || tau>=1
        error('Input tau must be a real positive scalar value.');
    end
end

end