function [D, x] = cheb_diff(N)        
   % create chebyshev points
   if N == 0
       x = 1;
       D = 0;
   end
   cheb_x = zeros(1, (N+1));
   for i = 1:N+1
       cheb_x(i) = cos((i-1)*pi/N);
   end
   
   
   % make differentiation matrix
   D = zeros(N+1, N+1);
   for i = 0:N
       for j = 0:N
           if i == j
               if i == 0 | i == N
                   D(1, 1) = (2*(N^2)+ 1)/(6);
                   D(N+1, N+1) = -1*D(1,1);
               else
                   D(i+1, j+1) = -1*cheb_x(j+1)/(2*(1 - (cheb_x(j+1))^2));
               end
           else
               if i == 0 | i == N c_i = 2;, else c_i = 1; end
               if j == 0 | j == N c_j = 2;, else c_j = 1; end
               D(i+1, j+1) = (c_i/c_j)*((-1)^((i+1)+(j+1)))/(cheb_x(i+1) - cheb_x(j+1));
           end
       end
   end
   
   x = cheb_x;
           
end