classdef aircraft
    properties
        m = 4.5;
        rho = 1.225;
        S = 0.473;
        g = 9.806;
        Cd0 = 0.0173;
        Cd1 = -0.0337;
        Cd2 = 0.0517;
        b = 3;
        p = 1;
        % limits is of the form [Clmin, Clmax, Vmin, Vmax, nu_min, nu_max, CTmin, CTmax, gamma_min, gamma_max]
        limits = [-0.2, 1.1, 10, 80, -pi/3, pi/3, -1e-4, 1e-4, -pi/4, pi/4]
        x
        u
        traj_params = struct();
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
            z = sigma(3);
            zdot = sigma(6); xdot = sigma(4); ydot = sigma(5); 
            zddot = sigma(9); xddot = sigma(7); yddot = sigma(8);
            
            % wind model
            p_exp = obj.p;
            Wx = VR*(-z)^p_exp;
            Wxz = (p_exp*VR)*((-z)^p_exp)/z;
            
            % non flat outputs
            V = ((xdot - Wx)^2 + ydot^2 + zdot^2)^0.5;
            Vdot = (xdot*xddot - xdot*zdot*Wxz - xddot*Wx + Wx*Wxz*zdot + ydot*yddot + zdot*zddot)/V;
            gamma = -asin(zdot/V);
            gammadot = (zdot*Vdot - V*zddot)/(V*(V^2 - zdot^2)^0.5);
            chi = atan2(ydot,(xdot - Wx));
            chidot = (xdot*yddot - yddot*Wx - ydot*xddot + ydot*zdot*Wxz)/(ydot^2 + xdot^2 + Wx^2 - 2*xdot*Wx);
            nu = atan((V*cos(gamma)*chidot - Wxz*zdot*sin(chi))/(V*gammadot + obj.g*cos(gamma) - Wxz*cos(chi)*sin(gamma)*zdot));
            Cl = (obj.m*V*cos(gamma)*chidot - obj.m*Wxz*zdot*sin(chi))/(0.5*obj.rho*obj.S*sin(nu)*V^2);
            if(nu == 0 && (chidot == 0 || zdot == 0 || chi == 0))
                Cl = 0;
            end
            
            % aerodynamic forces
            Cd = obj.Cd0 + obj.Cd1*Cl + obj.Cd2*Cl^2;
            D = 0.5*obj.rho*obj.S*V^2*Cd;
            T = obj.m*Vdot + D + obj.m*obj.g*sin(gamma) + obj.m*Wxz*zdot*cos(gamma)*cos(chi);
            CT = T/(0.5*obj.rho*obj.S*V^2);
            obj.x = [V, gamma, chi]; obj.u = [Cl, nu, CT];
        end
        
        function ydot = model(obj, X, U, VR)
            V = X(1,1); gamma = X(3,1); chi = X(2,1);
            z = X(6,1); 
            Cl = U(1,1); nu = U(2,1); CT = U(3,1); 
            
            Cd = obj.Cd0 + obj.Cd1*Cl + obj.Cd2*Cl^2;
            D = 0.5*obj.rho*obj.S*V^2*Cd;
            L = 0.5*obj.rho*obj.S*V^2*Cl;
            T = 0.5*obj.rho*obj.S*V^2*CT;
            
             % wind model
            p_exp = obj.p;
            Wx = VR*(-z)^p_exp;
            Wxz = (p_exp*VR)*((-z)^p_exp)/z;
            
            zdot = -V*sin(gamma); ydot = V*sin(chi)*cos(gamma); xdot = V*cos(chi)*cos(gamma) + Wx;
            Vdot = (-D/obj.m) - obj.g*sin(gamma) - Wxz*zdot*cos(chi)*cos(gamma) + (T/obj.m);
            chidot = (L*sin(nu)/(obj.m*V*cos(gamma))) + Wxz*zdot*sin(chi)/(V*cos(gamma));
            gammadot = (L*cos(nu)/(obj.m*V)) - obj.g*cos(gamma)/V + Wxz*zdot*cos(chi)*sin(gamma)/V;
            
            ydot = [Vdot; chidot; gammadot; xdot; ydot; zdot];
        end
        
        function A = get_jac(obj, t, tf, VR, coeffs, N)
            
            sigma = obj.get_traj(t, tf, coeffs, N);
            obj = obj.get_xu(sigma, VR);
            V = obj.x(1); gamma = obj.x(2); chi = obj.x(3);
            Cl = obj.u(1); nu = obj.u(2); CT = obj.u(3);
            
            % wind model
            p_exp = obj.p; z = sigma(3); zdot = sigma(6);
            %Wx = VR*(-z)^p_exp;
            Wxz = (p_exp*VR)*((-z)^p_exp)/z;
            
            
            % forces and drag polar
            Cd = obj.Cd0 + obj.Cd1*Cl + obj.Cd2*Cl^2;
            L = 0.5*obj.S*obj.rho*Cl*V^2;
            
            % constructing the Jacobian
            f11 = (obj.rho*obj.S*V*Cd/obj.m) - (obj.rho*obj.S*V*CT/obj.m);
            f12 = Wxz*zdot*sin(chi)*cos(gamma); f13 = -obj.g*cos(gamma) + Wxz*zdot*cos(chi)*sin(gamma);
            f14 = 0; f15 = 0; f16 = 0;
            f21 = (sin(nu)/(obj.m*cos(gamma)))*(0.5*obj.rho*obj.S*Cl) - (Wxz*zdot*sin(chi)/cos(gamma))*(1/V^2);
            f22 = (Wxz*zdot/(V*cos(gamma)))*cos(chi);
            f23 = (L*sin(nu)/(obj.m*V))*(tan(gamma)*sec(gamma)) + (Wxz*zdot*sin(chi)/V)*(tan(gamma)*sec(gamma));
            f24 = 0; f25 = 0; f26 = 0;
            f31 = (cos(nu)/obj.m)*(0.5*obj.rho*obj.S*Cl) + (obj.g*cos(gamma)/V^2) - (Wxz*zdot*cos(chi)*sin(gamma)/V^2);
            f32 = -Wxz*zdot*sin(chi)*cos(gamma)/V;
            f33 = (obj.g*sin(gamma) + Wxz*zdot*cos(gamma)*cos(chi))/V;
            f34 = 0; f35 = 0; f36 = 0;
            f41 = cos(chi)*cos(gamma); f42 = -V*cos(gamma)*sin(chi); f43 = -V*cos(chi)*sin(gamma); f44 = 0; f45 = 0; f46 = Wxz;
            f51 = sin(chi)*cos(gamma); f52 = V*cos(chi)*cos(gamma); f53 = -V*sin(gamma)*sin(chi);
            f54 = 0; f55 = 0; f56 = 0;
            f61 = -sin(gamma); f62 = 0; f63 = -V*cos(gamma);
            f64 = 0; f65 = 0; f66 = 0;
            
            A = [f11, f12, f13, f14, f15, f16;
                f21, f22, f23, f24, f25, f26;
                f31, f32, f33, f34, f35, f36;
                f41, f42, f43, f44, f45, f46;
                f51, f52, f53, f54, f55, f56;
                f61, f62, f63, f64, f65, f66];
         end
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
        