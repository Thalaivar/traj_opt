  function lam = Louiville(trajData, jacEval)
    intQuad = 0;
    for i = 1:trajData.N
        A = jacEval(trajData.X(i,:), trajData);
        intQuad = intQuad + trace(A);
    end
    
    lam = intQuad/trajData.N;
end