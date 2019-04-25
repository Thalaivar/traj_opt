function [D, x] = fourierDiff(N)
    h = 2*pi/N;
    x = linspace(h, 2*pi, N);
    col =  [0 0.5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]';
    D = toeplitz(col, col([1 N:-1:2]));
end