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