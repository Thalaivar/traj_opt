function FTM = get_FTM(aircraft, type)
    if strcmp(type, 'time_march')
        % FTM by explicit time evolution (N-pass)
         y_0 = eye(6); tspan = linspace(0, aircraft.tf,10000); y_T = zeros(6,6);
         for i = 1:6
             [~,y] = ode45(@(t,y) model(t, y, aircraft), tspan, y_0(:,i));
             y_T(:,i) = y(end,:)';
         end
         FTM = y_T;
    elseif strcmp(type, 'expo')
         % FTM by exponentials method
         t = linspace(aircraft.tf, 0, 1000);
         FTM = eye(6); del_t = t(1)-t(2);
         for i = 1:length(t)
            Jk = aircraft.get_jac(t(i), 'FD');
            FTM = FTM*expm(Jk*del_t);
         end
    elseif strcmp(type, 'friedmann')
        t = linspace(0, aircraft.tf, 1000);
        h = t(2) - t(1);
        FTM = eye(6);
        for i = 1:length(t)-2
            K = friedmann_K(h, (t(end) - i*h), aircraft);
            FTM = FTM*K;
        end
    end
     % model linearised about nominal trajectory
     function ydot = model(t, y, aircraft)
        A = aircraft.get_jac(t, 'FD');
        ydot = A*y;
     end    
end

