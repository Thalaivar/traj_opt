function visualisation(choice, aircraft, N)
    tf = aircraft.traj_params.tf;
    coeffs = aircraft.traj_params.coeffs;
    VR = aircraft.traj_params.VR;
    t = linspace(0,tf,1000);
    if choice == "state"
        X = zeros(length(t),3);
        for i = 1:length(t)
            sigma = aircraft.get_traj(t(i), tf, coeffs, N);
            aircraft = aircraft.get_xu(sigma, VR);
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
    end
    
    if choice == "traj-separate"
        X = zeros(length(t),3);
        for i = 1:length(t)
            sigma = aircraft.get_traj(t(i), tf, coeffs, N);
            X(i,:) = [sigma(1),sigma(2),sigma(3)];
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
    end
    
    if choice == "traj-3d"
        X = zeros(length(t),3);
        for i = 1:length(t)
            sigma = aircraft.get_traj(t(i), tf, coeffs, N);
            X(i,:) = [sigma(1),sigma(2),sigma(3)];
        end
        
        plot3(X(:,1), X(:,2), X(:,3))
        hold on
        plot3(X(1,1), X(1,2), X(1,3), 'sm')
        plot3(X(length(t)/4,1), X(length(t)/4,2), X(length(t)/4,3), 'or');
        plot3(X(floor(length(t)/8),1), X(floor(length(t)/8),2), X(floor(length(t)/8),3), 'og');
        xlabel('x'); ylabel('y'); zlabel('z');
        grid on
        axis equal
    end
end