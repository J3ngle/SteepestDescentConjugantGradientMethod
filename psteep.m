function [x,r,k,p,r_new,z] = psteep(A, P, b, itermax, tol)
    % Initialize variables
    x = zeros(size(b));
    r = b - A * x;
    z = P \ r;
    z=backwardSubstitution(P, r);
    p = z;
    
    
    for k = 1:itermax
        Ap = A * p;
        alpha = (z' * r) / (p' * Ap);
        x = x + alpha * p;
        r_new = r - alpha * Ap;
        
        residuals{k}=r_new;
        
        if norm(r_new) < tol
            break;
        end
        
        z_new = P \ r_new;
        beta = (z_new' * r_new) / (z' * r);
        p = z_new + beta * p;
        r = r_new;
        z = z_new;
    end
    
end
