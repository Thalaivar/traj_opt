function [FTM_expo, FTM_timemarch] = get_FTM(aircraft)
     % FTM by explicit time evolution (N-pass)
     y_0 = eye(6); tspan = linspace(0, aircraft.tf,10000); y_T = zeros(6,6);
     for i = 1:6
         [~,y] = ode45(@(t,y) model(t, y, aircraft), tspan, y_0(:,i));
         y_T(:,i) = y(end,:)';
     end
     FTM_timemarch = y_T;
     % FTM by exponentials method
     t = linspace(aircraft.tf, 0, 1000);
     FTM_expo = eye(6); del_t = t(1)-t(2);
     for i = 1:length(t)
        Jk = aircraft.get_jac(t(i), 'FD');
        FTM_expo = FTM_expo*expm(Jk*del_t);
     end
     % model linearised about nominal trajectory
     function ydot = model(t, y, aircraft)
        A = aircraft.get_jac(t, 'FD');
        ydot = A*y;
     end    
end
