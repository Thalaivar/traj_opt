function xguess = get_init_guess(type, N)
    if strcmp(type, 'circle')
        VR_0 = 0.1; tf_0 = 10;
        t = linspace(0, tf_0, 1000);
        x = 20*(-1 + cos(2*pi*t/tf_0));
        y = -20*(2^0.5)*sin(2*pi*t/tf_0);
        z = -20*(1 - cos(2*pi*t/tf_0)) - 0.11;    
    elseif strcmp(type, 'eight')
        VR_0 = 0.3; tf_0 = 20;
        t = linspace(0, tf_0, 1000);
        x = 20*sin(4*pi*t/tf_0);                  
        y = 40*sqrt(2)*sin(2*pi*t/tf_0);
        z = 20*(sin(4*pi*t/tf_0)+1) + 0.11; 
    else
        "ERROR: Incorrect choice for initial guess!!! Check in optimize_traj"
    end
    
    coeffs = get_coeffs([x',y',z'], tf_0, N);
    xguess = [coeffs(:,1); coeffs(:,2); coeffs(:,3)];
    xguess(3*(2*N+1)+1,1) = VR_0; xguess(3*(2*N+1)+2,1) = tf_0;
end