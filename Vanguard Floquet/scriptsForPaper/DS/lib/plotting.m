function plotting(X, p, type)
    if strcmp(type, 'circle')
        N = (length(X)-3)/8;
    else
        N = (length(X)-2)/8;
    end
    [~,x] = fourierdiff(N);
    t = X(8*N+1,1)*x/(2*pi);
    
    if strcmp(type, 'circle')
        chiLinearTerm = X(8*N+3);
        for i = 1:N
            X(N+i) = X(N+i) + chiLinearTerm*t(i);
        end
    end
    if p == 1
        profile = 'linear';
    else
        profile = 'exponential';
    end
    subplot(311)
    plot(t, X(1:N,1)); xlabel('t (s)'); ylabel('$V$', 'Interpreter', 'latex');
    grid minor; title(['State time history for N = ', num2str(N), ' in ', profile, ' wind'])
    subplot(312)
    plot(t, X(N+1:2*N,1)); xlabel('t (s)'); ylabel('$\chi$', 'Interpreter', 'latex');
    grid minor;
    subplot(313)
    plot(t, X(2*N+1:3*N,1)); xlabel('t (s)'); ylabel('$\gamma$', 'Interpreter', 'latex');
    grid minor;
    figure
    plot3(X(3*N+1:4*N,1), -X(4*N+1:5*N,1), -X(5*N+1:6*N,1), '--om'); xlabel('x'); ylabel('y'); zlabel('z');
    title(['Trajectory for N = ', num2str(N), ' in ', profile, ' wind'])
    grid minor;
    figure
    subplot(211)
    plot(t, X(6*N+1:7*N,1)); xlabel('t (s)'); ylabel('$C_L$', 'Interpreter', 'latex');
    grid minor; title(['Control time history for N = ', num2str(N), ' in ', profile, ' wind'])
    subplot(212)
    plot(t, X(7*N+1:8*N,1)); xlabel('t (s)'); ylabel('$\mu$', 'Interpreter', 'latex');
    grid minor;
end