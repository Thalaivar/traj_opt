function f = cheb_interp(cheb_x, x_k, N)
cheb_x
    poly_params = polyfit(cheb_x, x_k, N);
    x = linspace(cheb_x(1), cheb_x(N+1), (cheb_x(N+1) - cheb_x(1))*1000);
    f = zeros('like', x);
    for i = 1:length(x)
        temp = 0;
        for j = 1:length(poly_params) 
           temp = temp + poly_params(j)*(x(i)^(N+1 - j));
        end
        f(i) = temp;
    end
end