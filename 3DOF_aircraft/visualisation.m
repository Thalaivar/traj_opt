function visualisation(choice, aircraft)
    addpath('trajectory');
    
    tf = aircraft.tf; coeffs = aircraft.coeffs; N = aircraft.N;
    t = linspace(0,aircraft.tf,1000);
    if strcmp(choice, 'state')
        X = zeros(length(t),3);
        for i = 1:length(t)
            sigma = get_traj(t(i), tf, coeffs, N);
            aircraft = aircraft.get_xu(sigma);
            X(i,:) = aircraft.x;
        end
        
        subplot(3,1,1)
        plot(t,X(:,1))
        xlabel('t'); ylabel('V');
        grid on
        
        subplot(3,1,2)
        plot(t,X(:,3))
        xlabel('t'); ylabel('chi');
        grid on
        
        subplot(3,1,3)
        plot(t,X(:,2))
        xlabel('t'); ylabel('gamma');
        grid on
    
    elseif strcmp(choice,'traj-separate')
        X = zeros(length(t),3);
        for i = 1:length(t)
            sigma = get_traj(t(i), tf, coeffs, N);
            X(i,:) = [sigma(1),-sigma(2),-sigma(3)];
        end
        
        subplot(3,1,1)
        plot(t, X(:,1))
        xlabel('t'); ylabel('x');
        grid on
        
        subplot(3,1,2)
        plot(t, X(:,2))
        xlabel('t'); ylabel('y');
        grid on
        
        subplot(3,1,3)
        plot(t, X(:,3))
        xlabel('t'); ylabel('z');
        grid on
        
    elseif strcmp(choice, 'traj-3d')
        X = zeros(length(t),3);
        for i = 1:length(t)
            sigma = get_traj(t(i), tf, coeffs, N);
            X(i,:) = [sigma(1),-sigma(2),-sigma(3)];
        end
        
        comet3(X(:,1), X(:,2), X(:,3))
        xlabel('x'); ylabel('y'); zlabel('z');
        grid on 
        axis equal
    end
    
    rmpath('trajectory');
end