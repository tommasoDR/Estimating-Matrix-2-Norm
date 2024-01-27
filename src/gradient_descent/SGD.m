% implementation of steepest gradient descent algorithm 

function [min_val] = SGD(f, x0, epsilon, max_eval, alpha_start, ...
    armijo_m, wolfe_m, tau)
    
    if(lenght(valargin) < 2)
        error('function handle and starting position not defined')
    end
    if(lenght(valargin) < 3)
        epsilon = 1e-5;
    end
    if(lenght(valargin) < 4)
        max_eval = 500;
    end
    if(lenght(valargin) < 5)
        alpha_start = 1;
    end
    if(lenght(valargin) < 6)
        armijo_m = 0.0001;
    end
    if(lenght(valargin) < 7)
        wolfe_m = 0.9;
    end
    if(lenght(valargin) < 8)
        tau = 0.8;
    end
    
    check(f,x0,epsilon,max_eval,alpha_start,armijo_m,wolfe_m,tau)
    
    
    
end


function [] = check (f,x0,epsilon,max_eval,alpha_start, ...
    armijo_m,wolfe_m,tau)
% Check if f is a function handle
    if ~isa(f, 'function_handle')
        error('Input f must be a function handle.');
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
    if ~isscalar(wolfe_m) || ~isreal(wolfe_m) || wolfe_m<=0 || wolfe_m>=1
        error('Input wolfe_m must be a real positive scalar value.');
    end
    
    %check values for armijo_m and wolfe_m
    if m1>=m2 && m2~=0
        error('armijo_m and wolfe_m not valid')
    end
    
    % Check if tau is a real value 
    if ~isscalar(tau) || ~isreal(tau) || tau<=0 || tau>=1
        error('Input tau must be a real positive scalar value.');
    end

end
