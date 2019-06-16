function Xnew = expander(Xold,Nnew)
    Nold = (length(Xold)-3)/16;
    Xnew(16*Nnew+3,1) = 0;
    gridnew = linspace(0,1,Nnew); gridold = linspace(0,1,Nold);
    for i = 1:16
        Xnew((i-1)*Nnew+1:i*Nnew) = interp1(gridold,Xold((i-1)*Nold+1:i*Nold),gridnew);
    end
    Xnew(end-2:end) = Xold(end-2:end);

end