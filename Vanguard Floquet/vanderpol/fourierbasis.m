function coeffs = fourierbasis(X, T, N)
    % X has colums as time-series
    dim = size(X);
    M = dim(1); d = dim(2);% no. of timepoints
    t = linspace(0, T, M);
    phi = ones(M, 2*N+1);
    
    % least squares fit to find coeffs
    for i = 1:M
        for j = 2:N+1
            phi(i,j) = cos(2*pi*(j-1)*t(i)/T);
            phi(i,j+N) = sin(2*pi*(j-1)*t(i)/T);
        end
    end
    coeffs = zeros(2*N+1,d);
    for i = 1:d
        coeffs(:,i) = (inv(phi'*phi))*(phi'*X(:,i));
    end
end