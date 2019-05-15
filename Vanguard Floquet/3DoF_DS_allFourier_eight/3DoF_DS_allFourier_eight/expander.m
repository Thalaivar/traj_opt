function Xnew = expander(Xold,Nnew)
    Nold = (length(Xold)-2)/8;
    Xnew(8*Nnew+2,1) = 0;
    gridnew = linspace(0,1,Nnew); gridold = linspace(0,1,Nold);
    for i = 1:8
        Xnew((i-1)*Nnew+1:i*Nnew) = interp1(gridold,Xold((i-1)*Nold+1:i*Nold),gridnew);
    end
    Xnew(end-1:end) = Xold(end-1:end);

end