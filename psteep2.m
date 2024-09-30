function [x, residuals] = psteep2(A, P, b, intermax, tol)
    % Initialize variables
    x = zeros(size(b));
    r = b - A * x;
    z = P \ r;
    p = z;
    
    residuals = cell(intermax, 1); % Store residuals in a cell array
    
    for k = 1:intermax
        Ap = A * p;
        alpha = (z' * r) / (p' * Ap);
        x = x + alpha * p;
        r_new = r - alpha * Ap;
        
        residuals{k} = r_new; % Store the current residual
        
        if norm(r_new) < tol
            r = r_new; % Update r if we've converged
            break;
        end
        
        z_new = P \ r_new;
        beta = (z_new' * r_new) / (z' * r);
        p = z_new + beta * p;
        r = r_new;
        z = z_new;
    end
    
    % Trim the residuals cell array to the actual number of iterations
    residuals = residuals(1:k);
end
