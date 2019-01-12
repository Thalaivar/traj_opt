function FTM = get_FTM(solution, type)
    if strcmp(type, 'time_march')
        % FTM by explicit time evolution (N-pass)
         y_0 = eye(6); tspan = linspace(0, solution.tf,10000); y_T = zeros(6,6);
         for i = 1:6
             [~,y] = ode45(@(t,y) model(t, y, solution), tspan, y_0(:,i));
             y_T(:,i) = y(end,:)';
         end
         FTM = y_T;
    elseif strcmp(type, 'expo')
         % FTM by exponentials method
         t = linspace(solution.tf, 0, 1000);
         FTM = eye(3); del_t = t(1)-t(2);
         for i = 1:length(t)
            Jk = get_jac(solution, t(i), 'FD'); Jk = Jk(1:3,1:3);
            FTM = FTM*expm(Jk*del_t);
         end
    end
     % model linearised about nominal trajectory
     function ydot = model(t, y, solution)
        A = get_jac(solution, t, 'FD');
        ydot = A*y;
     end    
end

