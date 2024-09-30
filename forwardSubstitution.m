%forward substitution
function xtild = forwardSubstitution(L, b)
    % Check if L is a square matrix
    [m, n] = size(L);
    if m ~= n
        error('Matrix L must be square.');
    end
    
    % Check if L is lower triangular
    if ~istril(L)
        error('Matrix L must be lower triangular.');
    end
    
    % Initialize the solution vector
    xtild = zeros(n, 1);
    
    % Perform forward substitution
    for i = 1:n
        xtild(i) = (b(i) - L(i, 1:i-1) * xtild(1:i-1)) / L(i, i);
    end
end
