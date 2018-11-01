function SCmat = fBasis(t,tf,N)
% 26/12/17 12:18
% up to second derivative
% up to nth harmonic
% SCmat: 3x(2N+1)
temp(3,2*N) = 0;

    for i = 2*(1:N)
        mulplier = (pi*i/tf);
        temp(1,i-1) =  cos(mulplier*t);        temp(1,i) = sin(mulplier*t);
        temp(2,i-1) = -mulplier*temp(1,i);     temp(2,i) = mulplier*temp(1,i-1);
        temp(3,i-1) = -mulplier*temp(2,i);     temp(3,i) = mulplier*temp(2,i-1);
    end
    
SCmat = horzcat([1;0;0],temp);
end