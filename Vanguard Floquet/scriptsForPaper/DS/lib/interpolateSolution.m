function sol = interpolateSolution(type, X, N, shape)
    if strcmp(type, 'full-state')
        if nargin ~= 4
            error('Shape of trajectory not specified for full state collocation!');
        else
            if strcmp(shape, 'eight')
                oN = (length(X)-2)/8;
                sol = zeros(8*N+2,1);
            elseif strcmp(shape, 'circle')
                oN = (length(X)-3)/8;
                sol = zeros(8*N+3,1); sol(8*N+3) = X(8*oN+3);
            end
        end
        nXU = 8;
    elseif strcmp(type, 'diff-flat')
        oN = (length(X)-2)/3; nXU = 3;
        sol = zeros(3*N+2,1);
    end
    
    T = X(nXU*oN+1); VR = X(nXU*oN+2);
    
    [~,ot] = fourierdiff(oN); ot = T*ot/(2*pi);
    [~,t] = fourierdiff(N); t = T*t/(2*pi);
    for i = 1:nXU
        j = (i-1)*N; k = (i-1)*oN;
        sol(j+1:j+N) = interp_sinc(ot, X(k+1:k+oN), t);
    end
    sol(nXU*N+1) = T; sol(nXU*N+2) = VR;
end