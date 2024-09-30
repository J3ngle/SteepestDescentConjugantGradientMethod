%Backward substitution Function
function xback = backwardSubstitution(U, y)
    % Check if U is a square matrix
    [m, n] = size(U);
    if m ~= n
        error('Matrix U must be square.');
    end
    % Initialize the solution vector
    xback = zeros(n, 1);
    
    % Perform backward substitution
    for i = n:-1:1
        xback(i) = (y(i) - U(i, i+1:end) * xback(i+1:end)) / U(i, i);
    end
end
