function [dx,ddx] = dftDerv(x,T)
   N = numel(x);
   fac = 2*pi/T;
   xhat = fft(x); dxhat =1i*[0:N/2-1 0 -N/2+1:-1]'.*xhat; % 1 is odd
   dx = fac*real(ifft(dxhat));
   if nargout>1
       ddxhat = ((1i*[0:N/2 -N/2+1:-1]').^2).*xhat; % 2 is even
       ddx = fac*fac*real(ifft(ddxhat));
   end
end