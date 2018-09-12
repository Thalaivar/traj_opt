function f = poly_interp(x, y, N, t)
    prm = polyfit(x, y, N);
    f = zeros(size(t));
    for i = 1:length(t)
        for j = 1:length(prm)
            f(i) = f(i) + prm(j)*(t(i)^(N+1 - j));
        end
    end
end