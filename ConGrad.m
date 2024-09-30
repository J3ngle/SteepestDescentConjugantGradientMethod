function [conx, coniter, conz_new,conr_new] = ConGrad(A, P, b, intermax, tol)
    % Initialize the solution vector and other variables
    global n 
    conx = zeros(n, 1);
    r = b - A * conx;
    z = P \ r;
    p = z;
    
    for iter = 1:intermax
        alpha = (r' * z) / (p' * A * p);
        conx = conx + alpha * p;
        conr_new = r - alpha * A * p;
        
        % Check for convergence
        if norm(conr_new) < tol
            break;
        end
        
        conz_new = P \ conr_new;
        beta = (conr_new' * conz_new) / (r' * z);
        p = conz_new + beta * p;
        
        r = conr_new;
        z = conz_new;
    end
    
    if iter == intermax
        warning('PCG did not converge within the specified number of iterations.');
    end
end
