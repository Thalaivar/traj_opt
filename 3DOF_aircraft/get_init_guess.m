function xguess = get_init_guess(type, N, tf, VR)
    xguess = zeros(3*(2*N+1)+2,1);
    t = linspace(0, tf, 1000);
    if type == 'circle'
        x = 20*(-1 + cos(2*pi*t/tf));
        y = -20*(2^0.5)*sin(2*pi*t/tf);
        z = -20*(1 - cos(2*pi*t/tf)) - 0.11;
    end
    
    coeffs = get_coeffs([x',y',z'], tf, N);
    xguess = [coeffs(:,1); coeffs(:,2); coeffs(:,3)];
    xguess(3*(2*N+1)+1,1) = VR; xguess(3*(2*N+1)+2,1) = tf;
end