function [t, y] = RK4(model, tspan, y0, h)
    %% arrays to hold solution
    t = zeros(1, floor((tspan(2) - tspan(1))/h) + 1);
    y = zeros(floor((tspan(2) - tspan(1))/h) + 1, length(y0));
    y(1,:) = y0'; t(1) = tspan(1);
    %% 4th order RK iterations
    curr_t = tspan(1); curr_y = y0;
    i = 2;
    while curr_t + h <= tspan(2)
        % for notation, refer to wiki page on Runge-Kutta methods        
        k1 = h*model(curr_t, curr_y);
        k2 = h*model(curr_t + 0.5*h, curr_y + 0.5*k1);
        k3 = h*model(curr_t + 0.5*h, curr_y + 0.5*k2);
        k4 = h*model(curr_t + h, curr_y + k3);
        
      
        curr_y = curr_y + (k1 + 2*k2 + 2*k3 + k4)/6;
        curr_t = curr_t + h;
        
        t(i) = curr_t;
        y(i,:) = curr_y'; 
        i = i + 1;
    end
end