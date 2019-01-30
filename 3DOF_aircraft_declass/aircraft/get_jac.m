%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to calculate jacobian                                 %
% INPUT :                                                        %
%       t        - time                                          %
%       type     - either 'analytic' or 'FD'                     %
%       solution - struct with fields 'VR', 'tf', 'coeffs', 'N'  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function A = get_jac(solution, t, type)
    % model constants
    m = 4.5; rho = 1.225; S = 0.473;
    g = 9.806; Cd0 = 0.0173; Cd1 = -0.0337;
    Cd2 = 0.0517;
    
    tf = solution.tf; VR = solution.VR;
    coeffs = solution.coeffs; N = solution.N;
    
    if strcmp(type, 'analytic') 
        sigma = get_traj(t, tf, coeffs, N);
        [X, u] = get_xu(sigma, VR);
        V = X(1); gamma = X(2); chi = X(3);
        Cl = u(1); nu = u(2); CT = u(3);

        % wind model
        p_exp = 1; z = sigma(3); zdot = sigma(6);
        %Wx = VR*(-z)^p_exp;
        Wxz = -1*(p_exp*VR)*((-z)^(p_exp-1));


        % forces and drag polar
        Cd = Cd0 + Cd1*Cl + Cd2*Cl^2;
        L = 0.5*S*rho*Cl*V^2;

        % constructing the Jacobian
        f11 = (CT - Cd)*rho*S*V/m;
        f12 = Wxz*zdot*sin(chi)*cos(gamma); f13 = -1*g*cos(gamma) + Wxz*zdot*cos(chi)*sin(gamma);
        f14 = 0; f15 = 0; f16 = p_exp*(p_exp-1)*VR*((-z)^(p_exp-2))*zdot*cos(chi)*cos(gamma);
        f21 = 0.5*rho*S*Cl*sin(nu)/(m*cos(gamma)) - (Wxz*zdot*sin(chi)/(cos(gamma)*V^2));
        f22 = Wxz*cos(chi)*zdot/(V*cos(gamma));
        f23 = 0.5*rho*S*Cl*V*sin(nu)*sin(gamma)/(m*(cos(gamma))^2) + Wxz*zdot*sin(chi)*sin(gamma)/(V*(cos(gamma))^2);
        f24 = 0; f25 = 0; f26 = p_exp*(p_exp-1)*VR*((-z)^(p_exp-2))*zdot*sin(chi)/(V*cos(gamma));
        f31 = 0.5*rho*S*Cl*cos(nu)/m + (g*cos(gamma)/V^2) - (Wxz*cos(chi)*sin(gamma)*zdot/V^2);
        f32 = -Wxz*zdot*sin(chi)*sin(gamma)/V;
        f33 = g*sin(gamma)/V + (Wxz*zdot*cos(gamma)*cos(chi)/V);
        f34 = 0; f35 = 0; f36 = p_exp*(p_exp-1)*VR*((-z)^(p_exp-2))*cos(chi)*sin(gamma)*zdot/V;
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
    elseif strcmp(type, 'FD')
        I = eye(6); A = zeros(6);
        for i = 1:6
            % get states at current nominal point
            sig = get_traj(t, tf, coeffs, N);
            [X, ~] = get_xu(sig, VR);
            h = 1e-4; % perturbation
            x0 = [X(1); X(3); X(2); sig(1); sig(2); sig(3)];
            f1 = non_flat_model(t, x0 - h*I(:,i), solution); f2 = non_flat_model(t, x0 + h*I(:,i), solution);
            A(:,i) = (f2-f1)/(2*h);
        end
    end
 end