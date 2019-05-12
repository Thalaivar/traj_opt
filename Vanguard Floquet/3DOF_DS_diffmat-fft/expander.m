function Xnew = expander(Xold,Nnew)
    Nold = (length(Xold)-2)/3;
    Xnew(3*Nnew+2,1) = 0;
    Xnew(1:Nold) = Xold(1:Nold);
    Xnew(Nnew+1:Nnew+Nold) = Xold(Nold+1:2*Nold);
    Xnew(2*Nnew+1:2*Nnew+Nold) = Xold(2*Nold+1:3*Nold);
    Xnew(end-1:end) = Xold(end-1:end);

end