function X = fourierseries(t, T, coeffs)
    d = size(coeffs);
    N = (d(1)-1)/2; d = d(2);
    X = zeros(2,d);
    X(1,:) = coeffs(1,:);
    for i = 2:N+1
        for j = 1:d
            X(1,j) = X(1,j) + coeffs(i,j)*cos(2*pi*(i-1)*t/T) + coeffs(i+N,j)*sin(2*pi*(i-1)*t/T);
            X(2,j) = X(2,j) - coeffs(i,j)*(2*pi*(i-1)/T)*sin(2*pi*(i-1)*t/T) + coeffs(i+N,j)*(2*pi*(i-1)/T)*cos(2*pi*(i-1)*t/T);
        end
    end
end