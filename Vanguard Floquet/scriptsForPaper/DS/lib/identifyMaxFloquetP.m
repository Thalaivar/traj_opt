function maxFE = identifyMaxFloquetP(eigE, N, T)
    if(0.25*N < 2*pi/T)
        error("N is not large enough, must be at least 8*pi/T")
    end
    maxFE = Inf;
    % truncate eigenvalues to chop off spurious modes at the top
    eigE = eigE((abs(imag(eigE)) <= 0.25*N*2*pi/T));
    % sort eigenvalues according to increasing real part
    sortedEigE = sort(eigE, 'descend', 'ComparisonMethod', 'real');
    lagVal = 5;
    for i = 1:length(sortedEigE)-lagVal
       if ~imModulo(sortedEigE(i:i+lagVal),2*pi/T)
           maxFE = real(sortedEigE(i));
           break
       end
    end
end
function prodz = imModulo(x,y)
    % test if the differences of elements of x are multiples of y:
    N = numel(x);
    if N<5
        error('N>=5 is required');
    end

    z = [];
    countr = 0;
    for i = 2:N
        z(N-i+1) = mod(imag(x(1))-imag(x(i)),y);
        if abs(z(N-i+1))<5e-4
           countr = countr+1;
        end
    end

    if countr>=2
       prodz = 0;
    else
       prodz = 1;
    end

end