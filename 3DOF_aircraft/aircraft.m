classdef aircraft
    properties
        m = 4.5;
        rho = 1.225;
        S = 0.473;
        g = 9.8;
        Cd0 = 0.0173;
        Cd1 = -0.0337;
        Cd2 = 0.0517;
        b = 3;
        p = 1;
        % limits is of the form [Clmin, Clmax, Vmin, Vmax, nu_min, nu_max, CTmin, CTmax, gamma_min, gamma_max]
        limits = [-0.2, 1.1, 10, 26, -pi/3, pi/3, -1e-4, 1e-4, -pi/4, pi/4]
        x
        u
    end
    
    methods
        % constructor
        function obj = aircraft(m,rho,S,g,Cd0,Cd1,Cd2,b,p,limits)
           if nargin > 0
                obj.m = m; obj.rho = rho; obj.S = S; obj.Cd0 = Cd0;
                obj.Cd1 = Cd1; obj.Cd2 = Cd2; obj.b = b; obj.p = p;
                obj.limits = limits;
           end 
        end
        
        % to get the non flat states
        function obj = get_xu(obj, sigma, VR)
            z = sigma(1);
            zdot = sigma(6); xdot = sigma(4); ydot = sigma(5); 
            zddot = sigma(9); xddot = sigma(7); yddot = sigma(8);
            
            % wind model
            p_exp = obj.p;
            Wx = VR*(-z)^p_exp;
            Wxz = (p_exp*VR)*((-z)^p_exp)/z;
            
            % non flat outputs
            V = ((xdot - Wx)^2 + ydot^2 + zdot^2)^0.5;
            Vdot = (xdot*xddot - xdot*zdot*Wxz - xddot*Wx + Wx*Wxz*zdot + ydot*yddot + zdot*zddot)/V;
            gamma = asin(-zdot/V);
            gammadot = (zdot*Vdot - V*zddot)/(V*(V^2 - zdot^2)^0.5);
            chi = atan(ydot/(xdot - Wx));
            chidot = (xdot*yddot - yddot*Wx - ydot*xddot + ydot*zdot*Wxz)/(ydot^2 + xdot^2 + Wx^2 - 2*xdot*Wx);
            nu = atan((V*cos(gamma)*chidot - Wxz*zdot*sin(chi))/(V*gammadot + obj.g*cos(gamma) - Wxz*cos(chi)*sin(gamma)*zdot));
            Cl = (obj.m*V*cos(gamma)*chidot - obj.m*Wxz*zdot*sin(chi))/(0.5*obj.rho*obj.S*sin(nu)*V^2);
            
            % aerodynamic forces
            Cd = obj.Cd0 + obj.Cd1*Cl + obj.Cd2*Cl^2;
            D = 0.5*obj.rho*obj.S*V^2*Cd;
            T = obj.m*Vdot + D + obj.m*obj.g*sin(gamma) + obj.m*Wxz*zdot*cos(gamma)*cos(chi);
            CT = T/(0.5*obj.rho*obj.S*V^2);
            obj.x = [V, gamma, chi]; obj.u = [Cl, nu, CT];
        end
s
%         function A = get_jac(obj, t, tf, a, b, VR)
%             sigma = get_traj(t, tf, a, b);
%             [x, u] = get_xu(obj, sigma, VR);
%             
%             % wind model
%             p_exp = obj.p;
%             Wx = VR*(-z)^p_exp;
%             Wxz = (p_exp*VR)*((-z)^p_exp)/z;
%             
%             % aerodynamic forces and drag polar
%             A = []
%         end
    end
    
    methods (Static)
        
        function sigma = get_traj(t, tf, coeffs, N)
            ax = coeffs(:,1); ay = coeffs(:,2); az = coeffs(:,3);
            
            x = ax(1); y = ay(1); z = az(1); 
            xdot = 0; ydot = 0; zdot = 0; xddot = 0; yddot = 0; zddot = 0; 
            for i = 2:N+1
                x = x + ax(i)*cos(2*pi*(i-1)*t/tf) + ax(i+N)*sin(2*pi*(i-1)*t/tf);
                y = y + ay(i)*cos(2*pi*(i-1)*t/tf) + ay(i+N)*sin(2*pi*(i-1)*t/tf);
                z = z + az(i)*cos(2*pi*(i-1)*t/tf) + az(i+N)*sin(2*pi*(i-1)*t/tf);
                xdot = xdot - ax(i)*(2*pi*(i-1)/tf)*sin(2*pi*(i-1)*t/tf) + ax(i+N)*(2*pi*(i-1)/tf)*cos(2*pi*(i-1)*t/tf);
                ydot = ydot - ay(i)*(2*pi*(i-1)/tf)*sin(2*pi*(i-1)*t/tf) + ay(i+N)*(2*pi*(i-1)/tf)*cos(2*pi*(i-1)*t/tf);
                zdot = zdot - az(i)*(2*pi*(i-1)/tf)*sin(2*pi*(i-1)*t/tf) + az(i+N)*(2*pi*(i-1)/tf)*cos(2*pi*(i-1)*t/tf);
                xddot = xddot - ax(i)*((2*pi*(i-1)/tf)^2)*cos(2*pi*(i-1)*t/tf) - ax(i+N)*((2*pi*(i-1)/tf)^2)*sin(2*pi*(i-1)*t/tf);
                yddot = yddot - ay(i)*((2*pi*(i-1)/tf)^2)*cos(2*pi*(i-1)*t/tf) - ay(i+N)*((2*pi*(i-1)/tf)^2)*sin(2*pi*(i-1)*t/tf);
                zddot = zddot - az(i)*((2*pi*(i-1)/tf)^2)*cos(2*pi*(i-1)*t/tf) - az(i+N)*((2*pi*(i-1)/tf)^2)*sin(2*pi*(i-1)*t/tf);
            end
            
            sigma = [x,y,z,xdot,ydot,zdot,xddot,yddot,zddot];
        end
        
    end
end      
        