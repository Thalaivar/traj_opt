function [ixdot, iydot] = integrate_xy(state, tf, VR, p)
    N = size(state); N = N(1);
    ixdot = 0; iydot = 0;
    [~,x] = fourierdiff(N);
    V = state(:,1);
    chi = state(:,2);
    gamma = state(:,3);
    z = state(:,4);
    t = tf*x/(2*pi);
    h = t(2)-t(1);
    for i = 1:N-1
        Wx1 = VR*(-z(i))^p; Wx2 = VR*(-z(i+1))^p;
        f_ixdot = V(i)*cos(chi(i))*cos(gamma(i)) + Wx1 + V(i+1)*cos(chi(i+1))*cos(gamma(i+1)) + Wx2;
        ixdot = ixdot + 0.5*h*f_ixdot;
        f_iydot = V(i)*sin(chi(i))*cos(gamma(i)) + V(i+1)*sin(chi(i+1))*cos(gamma(i+1));
        iydot = iydot + 0.5*h*f_iydot;
    end
end