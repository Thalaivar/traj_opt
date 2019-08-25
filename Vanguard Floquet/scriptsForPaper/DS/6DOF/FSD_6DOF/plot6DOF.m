function plot6DOF(sol, N, p, shape)
    nState = 12; nControl = 4; nOptim = nState + nControl;
    % get states and controls
    X = zeros(N,nState); U = zeros(N,nControl);
    for i = 1:nOptim
        j = (i-1)*N;
        if i <= nState
            X(:,i) = sol(j+1:j+N);
        else
            U(:,i-nState) = sol(j+1:j+N);
        end
    end
    
    T = sol(nOptim*N+1); VR = sol(nOptim*N+2);
    [~,fourierGrid] = fourierdiff(N);
    t = T*fourierGrid/(2*pi);
    
    % add linear term to psi for loiter
    if strcmp(shape, 'loiter')
        psiLinearterm = sol(nOptim*N+3);
        for i = 1:N
            X(i,9) = X(i,9) + psiLinearterm*t(i);
        end
    end
    
    if p < 1, windType = 'exponential'; else, windType = 'linear'; end
    
    % plot trajectory
    figure
    plot3(X(:,10), -X(:,11), -X(:,nState), '--om')
    title([shape,' trajectory in ',windType,' profile; N = ',num2str(N)]);
    xlabel('x'); ylabel('y'); zlabel('z');
    grid minor
    
    % plot euler angles
    figure
    subplot(311)
    plot(t, X(:,7))
    title(['Euler angles for ',shape,' trajectory in ',windType,' wind; N = ', num2str(N)]);
    xlabel('t(s)'); ylabel('$\phi$', 'Interpreter', 'latex');
    subplot(312)
    plot(t,X(:,8))
    xlabel('t(s)'); ylabel('$\theta$', 'Interpreter', 'latex');
    subplot(313)
    plot(t,X(:,9))
    xlabel('t(s)'); ylabel('$\psi$', 'Interpreter', 'latex');
    
    % plot velocities
    figure
    subplot(311)
    plot(t, X(:,1))
    title(['Velocities for ',shape,' trajectory in ',windType,' wind; N = ', num2str(N)]);
    xlabel('t(s)'); ylabel('$u$', 'Interpreter', 'latex');
    subplot(312)
    plot(t,X(:,2))
    xlabel('t(s)'); ylabel('$v$', 'Interpreter', 'latex');
    subplot(313)
    plot(t,X(:,3))
    xlabel('t(s)'); ylabel('$w$', 'Interpreter', 'latex');
    
    % plot angular rates
    figure
    subplot(311)
    plot(t, X(:,4))
    title(['Angular rates for ',shape,' trajectory in ',windType,' wind; N = ', num2str(N)]);
    xlabel('t(s)'); ylabel('$p$', 'Interpreter', 'latex');
    subplot(312)
    plot(t,X(:,5))
    xlabel('t(s)'); ylabel('$q$', 'Interpreter', 'latex');
    subplot(313)
    plot(t,X(:,6))
    xlabel('t(s)'); ylabel('$r$', 'Interpreter', 'latex');
    
    % plot control inputs
    figure
    subplot(221)
    plot(t, U(:,1))
    xlabel('t(s)'); ylabel('$\delta_f$', 'Interpreter', 'latex');
    subplot(222)
    plot(t, U(:,2))
    xlabel('t(s)'); ylabel('$\delta_e$', 'Interpreter', 'latex');
    subplot(223)
    plot(t, U(:,3))
    xlabel('t(s)'); ylabel('$\delta_a$', 'Interpreter', 'latex');
    subplot(224)
    plot(t, U(:,4))
    xlabel('t(s)'); ylabel('$\delta_r$', 'Interpreter', 'latex');
end