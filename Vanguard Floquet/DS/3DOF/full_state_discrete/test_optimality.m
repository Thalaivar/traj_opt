function [c, ceq, X0] = test_optimality(ac, N)
    [~,x] = fourierdiff(N);
    t = ac.tf*x/(2*pi);
%     [~,x] = chebdiff(N-1);
%     t = 0.5*ac.tf*(1-x);
    state = zeros(N,6); control = zeros(N,3);
    for i = 1:N
        sig = get_traj(t(i), ac.tf, ac.coeffs, ac.N);
        state(i,4:6) = [sig(1), sig(2), sig(3)];
        ac = ac.get_xu(sig);
        V = ac.x(1); gamma = ac.x(2); chi = ac.x(3);
        Cl = ac.u(1); mu = ac.u(2); CT = ac.u(3);
        state(i,1:3) = [V, chi, gamma];
        control(i,:) = [Cl, mu, CT];
    end
    X0 = zeros(9*N+2,1);
    for i = 1:9  
        j = (i-1)*N;
        if i <= 6
            X0(j+1:j+N,1) = state(:,i);
        else
            X0(j+1:j+N,1) = control(:,i-6);
        end
    end
    X0(9*N+1,1) = ac.tf; X0(9*N+2,1) = ac.VR;
    
    [c, ceq] = constrain_function(X0, ac.p);
end