function plot_xyz(sol,p)
    N = (length(sol)-2)/7;
    V = sol(1:N,1); chi = sol(N+1:2*N,1); 
    gamma = sol(2*N+1:3*N,1); z = sol(3*N+1:4*N,1);
    T = sol(7*N+1,1); VR = sol(7*N+2,1);
    [~,x] = fourierdiff(N);
    t = T*x/(2*pi);
    h = t(2)-t(1);
    pos_x = zeros(N,1);
    pos_y = zeros(N,1);
    ixdot = 0; iydot = 0;
    for i = 1:N-1
        Wx1 = VR*(-z(i))^p; Wx2 = VR*(-z(i+1))^p;
        f_ixdot = V(i)*cos(chi(i))*cos(gamma(i)) + Wx1 + V(i+1)*cos(chi(i+1))*cos(gamma(i+1)) + Wx2;
        ixdot = ixdot + 0.5*h*f_ixdot;
        f_iydot = V(i)*sin(chi(i))*cos(gamma(i)) + V(i+1)*sin(chi(i+1))*cos(gamma(i+1));
        iydot = iydot + 0.5*h*f_iydot;
        pos_x(i+1) = ixdot; pos_y(i+1) = iydot;
    end
    plot3(pos_x, -pos_y, -z)
    xlabel('x'); ylabel('y'); zlabel('z');
    grid minor;
end