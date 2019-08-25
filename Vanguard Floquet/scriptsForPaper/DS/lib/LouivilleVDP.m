function lam = LouivilleVDP(trajData)
    intQuad = 0;
    for i = 1:trajData.N
        A = jacVDP(trajData.X(i,:), trajData);
        intQuad = intQuad + trace(A);
    end
    
    lam = intQuad/trajData.N;
end