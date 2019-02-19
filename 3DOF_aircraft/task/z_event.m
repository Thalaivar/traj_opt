function [z, isterminal, direction] = z_event(~,y)
    z = y(6,1)+1e-4;
    isterminal = 1;
    direction = 0;
end