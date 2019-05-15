classdef aircraft
%     3DOF aircraft object for 
%     used for 3DOF diff-flat diagnostics
    properties 
        m
        b
        S
        CD0
        CD1
        CD2
        % normalization
        NCT
        
        % wind exponent
        p_exp
    end
    methods
        % contructor
        function obj = aircraft(m,b,S,CD0,CD1,CD2,p_exp)
            obj.m = m;
            obj.b = b;
            obj.S = S;
            obj.CD0 = CD0;
            obj.CD1 = CD1;
            obj.CD2 = CD2;
            obj.p_exp = p_exp;
        end
        function dZ = stateDervs(obj,Z,VR)
            % Z = [V,chi,gam,x,y,z,CL,mu,CT]
            V = Z(1); chi = Z(2); gam = Z(3); % x =  Z(4); y = Z(5); 
            z = Z(6); CL = Z(7); mu = Z(8); CT = Z(9);
            
            dZ = zeros(6,1);
            
            %%%%% Wind Model
            % exponential
            Wx = VR*(-z)^obj.p_exp;
            Wxz = (obj.p_exp*VR)*((-z)^obj.p_exp)/z;  
            
            Cmu = cos(mu); Smu = sin(mu);
            Cgam = cos(gam); Sgam = sin(gam);
            Cchi = cos(chi); Schi = sin(chi);
            
            Q = 0.5*1.225*V*V;
            D = Q*obj.S*(obj.CD0 + obj.CD1*CL + obj.CD2*CL*CL);
            L = Q*obj.S*CL;
            T = Q*obj.S*CT;
            
            dZ(4) = V*Cchi*Cgam + Wx;
            dZ(5) = V*Schi*Cgam;
            dZ(6) = -V*Sgam;
            
            dZ(1) = -D/obj.m - 9.806*Sgam - Wxz*dZ(6)*Cchi*Cgam + T/obj.m;
            dZ(2) = L*Smu/(obj.m*V*Cgam) + Wxz*dZ(6)*Schi/(V*Cgam);
            dZ(3) = L*Cmu/(obj.m*V) - 9.806*Cgam/V + Wxz*Cchi*Sgam*dZ(6)/V;
        end    
    end
end