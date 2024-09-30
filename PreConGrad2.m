function [xcon, iter,z_new,p,rcon] = PreConGrad2(A, P, b, maxIter, tol)

    % Initialize the solution vector and other variables
    xcon = zeros(n, 1);
    rcon = b - A * xcon;
    z = P \ rcon;
    p = z;
    
    for iter = 1:maxIter
        alpha = (rcon' * z) / (p' * A * p);
        xcon = xcon + alpha * p;
        r_new = rcon - alpha * A * p;
        
        % Check for convergence
        if norm(r_new) < tol
            break;
        end
        
        z_new = P \ r_new;
        beta = (r_new' * z_new) / (rcon' * z);
        p = z_new + beta * p;
        
        rcon = r_new;
        z = z_new;
    end
    
    if iter == maxIter
        warning('PCG did not converge within the specified number of iterations.');
    end
end
