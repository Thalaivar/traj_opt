function fourCoeff = fitFS(Z,tf,M)
% determines the fourier coefficients when the period is known
% Z must be a tall matrix; time-series along the columns
    N = size(Z,1);
    t = linspace(0,tf,N)';
    A = ones(length(t),2*M+1);
    for i = 1:M
        A(:,[2*i 2*i+1]) = [cos((2*pi*t/tf)*i),sin((2*pi*t/tf)*i)];
    end
    fourCoeff = (A\Z)';
% fourCoeff has the coefficients arranged along the columns    
end