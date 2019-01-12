function coeffs = get_coeffs(X, tf, N)
    % X has colums as x,y,z time-series
    x = X(:,1); y = X(:,2); z = X(:,3);
    M = length(x); % no. of timepoints
    t = linspace(0, tf, M);
    phi = ones(M, 2*N+1);
    
    % least squares fit to find coeffs
    for i = 1:M
        for j = 2:N+1
            phi(i,j) = cos(2*pi*(j-1)*t(i)/tf);
            phi(i,j+N) = sin(2*pi*(j-1)*t(i)/tf);
        end
    end
    coeffs_x = (inv(phi'*phi))*(phi'*x);
    coeffs_y = (inv(phi'*phi))*(phi'*y);
    coeffs_z = (inv(phi'*phi))*(phi'*z);
    
    coeffs = [coeffs_x, coeffs_y, coeffs_z];
end