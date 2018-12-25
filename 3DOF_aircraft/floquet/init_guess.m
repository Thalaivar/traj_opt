% to generate initial guess for eigenvector collocation method of
% estimation of floquet exponents (probably obsolete)
% aircraft object needs to have tf, VR, coeffs and params set
function [init_guess, M] = init_guess(aircraft, FTM)
   M = 100; % no. of points in collocation grid (needs to be same during optimisation)
   [~, cheb_x] = cheb_diff(M-1);
   cheb_t = 0.5*aircraft.tf*(1 - cheb_x);

   [V,D] = eig(FTM);
   eigval = (1/aircraft.tf)*log(D(3,3)); 
   params = [real(eigval), imag(eigval)];
   y0 = V(:,3);
   tspan = cheb_t;

   [~, y] = ode45(@(t,y) model(t, y, params, aircraft), tspan, y0);

   init_guess = zeros(6*M,1);
   for i = 1:6
       j = (i-1)*M + 1;
       init_guess(j:j+M-1,1) = y(:,i);
   end

   init_guess(end+1,1) = real(eigval);
   init_guess(end+1,1) = imag(eigval);

   function ydot = model(t, y, params, aircraft)
        alpha = y(1:3,1); beta = y(4:6,1);
        eig_val = complex(params(1), params(2));
        A = aircraft.get_jac(t);
        I = eye(3);
        alpha_dot = (A - real(eig_val)*I)*alpha + imag(eig_val)*beta;
        beta_dot  = (A - real(eig_val)*I)*beta  - imag(eig_val)*alpha; 
        ydot = [alpha_dot;beta_dot];
   end
end