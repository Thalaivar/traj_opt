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
    
    % check error in result
    err = zeros(d,M);
    for ii = 1:M
        fit_result = zeros(1,3);
        fit_result(1,:) = coeffs(1,:); 
        for i = 2:N+1
            for j = 1:d
                fit_result(1,j) = fit_result(1,j) + coeffs(i,j)*cos(2*pi*(i-1)*t(ii)/T) + coeffs(i+N,j)*sin(2*pi*(i-1)*t(ii)/T);
            end
        end
        for j = 1:d
            err(j,ii) = (X(ii,j) - fit_result(1,j));
        end
    end
    rmse = 0;
    for i = 1:d
        rmse = rmse + norm(err(i,:))/sqrt(length(err(i,:)));
    end
    rmse = rmse/d;
    fprintf('The RMSE of the Fourier fit is: %0.4f\n', rmse);
end