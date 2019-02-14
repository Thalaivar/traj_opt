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
        x
        u
        coeffs
        tf
        VR
        N
        solution
    end
    
    methods
        % constructor
        function obj = aircraft(m,rho,S,Cd0,Cd1,Cd2,b,p)
           if nargin > 0
                obj.m = m; obj.rho = rho; obj.S = S; obj.Cd0 = Cd0;
                obj.Cd1 = Cd1; obj.Cd2 = Cd2; obj.b = b; obj.p = p;
           end 
        end
        
        % to get the non flat states
        function obj = get_xu(obj, sigma)
            z = sigma(3);
            zdot = sigma(6); xdot = sigma(4); ydot = sigma(5); 
            zddot = sigma(9); xddot = sigma(7); yddot = sigma(8);
            
            % wind model
            p_exp = obj.p;
            Wx = obj.VR*(-z)^p_exp;
            Wxz = (p_exp*obj.VR)*((-z)^p_exp)/z;
            
            % non flat outputs
            V = ((xdot - Wx)^2 + ydot^2 + zdot^2)^0.5;
            Vdot = (xdot*xddot - xdot*zdot*Wxz - xddot*Wx + Wx*Wxz*zdot + ydot*yddot + zdot*zddot)/V;
            gamma = -asin(zdot/V);
            gammadot = (zdot*Vdot - V*zddot)/(V*(V^2 - zdot^2)^0.5);
            chi = atan2(ydot,(xdot - Wx));
            chidot = (xdot*yddot - yddot*Wx - ydot*xddot + ydot*zdot*Wxz)/(ydot^2 + xdot^2 + Wx^2 - 2*xdot*Wx);
            nu = atan((V*cos(gamma)*chidot - Wxz*zdot*sin(chi))/(V*gammadot + obj.g*cos(gamma) - Wxz*cos(chi)*sin(gamma)*zdot));
            Cl = (obj.m*V*cos(gamma)*chidot - obj.m*Wxz*zdot*sin(chi))/(0.5*obj.rho*obj.S*sin(nu)*V^2);
%             if(nu == 0 && (chidot == 0 || zdot == 0 || chi == 0))
%                 Cl = 0;
%             end
            
            % aerodynamic forces
            Cd = obj.Cd0 + obj.Cd1*Cl + obj.Cd2*Cl^2;
            D = 0.5*obj.rho*obj.S*V^2*Cd;
            T = obj.m*Vdot + D + obj.m*obj.g*sin(gamma) + obj.m*Wxz*zdot*cos(gamma)*cos(chi);
            CT = T/(0.5*obj.rho*obj.S*V^2);
            obj.x = [V, gamma, chi]; obj.u = [Cl, nu, CT];
        end
        
        function ydot = non_flat_model(obj, t, X)
            V = X(1,1); gamma = X(3,1); chi = X(2,1);
            z = X(6,1);  
                
            sigma = get_traj(t, obj.tf, obj.coeffs, obj.N);
            obj = obj.get_xu(sigma);
            Cl = obj.u(1); nu = obj.u(2); CT = obj.u(3);
            
            Cd = obj.Cd0 + obj.Cd1*Cl + obj.Cd2*Cl^2;
            D = 0.5*obj.rho*obj.S*V^2*Cd;
            L = 0.5*obj.rho*obj.S*V^2*Cl;
            T = 0.5*obj.rho*obj.S*V^2*CT;
            
             % wind model
            p_exp = obj.p;
            Wx = obj.VR*(-z)^p_exp;
            Wxz = (p_exp*obj.VR)*((-z)^p_exp)/z;
            
            zdot = -V*sin(gamma); xdot = V*cos(chi)*cos(gamma) + Wx; ydot = V*sin(chi)*cos(gamma);
            Vdot = (-D/obj.m) - obj.g*sin(gamma) - Wxz*zdot*cos(chi)*cos(gamma) + (T/obj.m);
            chidot = (L*sin(nu)/(obj.m*V*cos(gamma))) + Wxz*zdot*sin(chi)/(V*cos(gamma));
            gammadot = (L*cos(nu)/(obj.m*V)) - obj.g*cos(gamma)/V + Wxz*zdot*cos(chi)*sin(gamma)/V;
            
            ydot = [Vdot; chidot; gammadot; xdot; ydot; zdot];
            if(imag(ydot) ~= 0)
                ydot
            end
        end
        
        function A = get_jac(obj, t, type)
            if strcmp(type, 'analytic')
               % if obj.p == 1  
                    sigma = get_traj(t, obj.tf, obj.coeffs, obj.N);
                    obj = obj.get_xu(sigma);
                    V = obj.x(1); gamma = obj.x(2); chi = obj.x(3);
                    Cl = obj.u(1); nu = obj.u(2); CT = obj.u(3);

                    % wind model
                    p_exp = obj.p; z = sigma(3); zdot = sigma(6);
                    %Wx = VR*(-z)^p_exp;
                    Wxz = -1*(p_exp*obj.VR)*((-z)^(p_exp-1));


                    % forces and drag polar
                    Cd = obj.Cd0 + obj.Cd1*Cl + obj.Cd2*Cl^2;
                    L = 0.5*obj.S*obj.rho*Cl*V^2;

                    % constructing the Jacobian
                    f11 = (CT - Cd)*obj.rho*obj.S*V/obj.m;
                    f12 = Wxz*zdot*sin(chi)*cos(gamma); f13 = -1*obj.g*cos(gamma) + Wxz*zdot*cos(chi)*sin(gamma);
                    f14 = 0; f15 = 0; f16 = p_exp*(p_exp-1)*obj.VR*((-z)^(p_exp-2))*zdot*cos(chi)*cos(gamma);
                    f21 = 0.5*obj.rho*obj.S*Cl*sin(nu)/(obj.m*cos(gamma)) - (Wxz*zdot*sin(chi)/(cos(gamma)*V^2));
                    f22 = Wxz*cos(chi)*zdot/(V*cos(gamma));
                    f23 = 0.5*obj.rho*obj.S*Cl*V*sin(nu)*sin(gamma)/(obj.m*(cos(gamma))^2) + Wxz*zdot*sin(chi)*sin(gamma)/(V*(cos(gamma))^2);
                    f24 = 0; f25 = 0; f26 = p_exp*(p_exp-1)*obj.VR*((-z)^(p_exp-2))*zdot*sin(chi)/(V*cos(gamma));
                    f31 = 0.5*obj.rho*obj.S*Cl*cos(nu)/obj.m + (obj.g*cos(gamma)/V^2) - (Wxz*cos(chi)*sin(gamma)*zdot/V^2);
                    f32 = -Wxz*zdot*sin(chi)*sin(gamma)/V;
                    f33 = obj.g*sin(gamma)/V + (Wxz*zdot*cos(gamma)*cos(chi)/V);
                    f34 = 0; f35 = 0; f36 = p_exp*(p_exp-1)*obj.VR*((-z)^(p_exp-2))*cos(chi)*sin(gamma)*zdot/V;
                    f41 = cos(chi)*cos(gamma); f42 = -V*sin(chi)*cos(gamma); f43 = -V*cos(chi)*sin(gamma);
                    f44 = 0; f45 = 0; f46 = Wxz;
                    f51 = sin(chi)*cos(gamma); f52 = V*cos(chi)*cos(gamma); f53 = -V*sin(chi)*sin(gamma);
                    f54 = 0; f55 = 0; f56 = 0;
                    f61 = -sin(gamma); f62 = 0; f63 = -V*cos(gamma); f64 = 0; f65 = 0; f66 = 0;
                    A = [f11, f12, f13, f14, f15, f16;
                         f21, f22, f23, f24, f25, f26;
                         f31, f32, f33, f34, f35, f36;
                         f41, f42, f43, f44, f45, f46;
                         f51, f52, f53, f54, f55, f56;
                         f61, f62, f63, f64, f65, f66];
                %else, "ERROR! Analytic Jacobian cannot be used for exponential profile!", A = [];
            elseif strcmp(type, 'FD')
                I = eye(6); A = zeros(6);
                for i = 1:6
                    % get states at current nominal point
                    sig = get_traj(t, obj.tf, obj.coeffs, obj.N);
                    obj = obj.get_xu(sig);
                    h = 1e-4; % perturbation
                    x0 = [obj.x(1); obj.x(3); obj.x(2); sig(1); sig(2); sig(3)];
                    f1 = obj.non_flat_model(t, x0 - h*I(:,i)); f2 = obj.non_flat_model(t, x0 + h*I(:,i));
                    A(:,i) = (f2-f1)/(2*h);
                end
            end
         end
    end
end      
        