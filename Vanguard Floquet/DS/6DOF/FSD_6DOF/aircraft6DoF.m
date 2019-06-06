classdef aircraft6DoF
 
   properties
       m % mass
       S %
       b % Reference quantities
       c %
       d %
       
       Ix  %
       Iy  % Moments of intertia about body axes
       Iz  %
       Ixz %
              
       % Coefficients and Derivatives
       CD0
       CD1
       CD2       
       CYb
       CYp
       CYr
       CYda
       CYdr
       CL0
       CLalf
       CLq
       CLdf
       CLde
       Clb
       Clp
       Clr
       Clda
       Cldr
       Cm0
       Cmalf
       Cmq
       Cmdf
       Cmde
       Cnb
       Cnp
       Cnr
       Cnda
       Cndr
       
       p_exp
   end
   methods
       % constructor
       function obj = aircraft6DoF(m,S,b,c,d,Ix,Iy,Iz,Ixz,coeff,p_exp)
           obj.m = m;
           obj.S = S; 
           obj.b = b; 
           obj.c = c;
           obj.d = d;
           
           obj.Ix = Ix;
           obj.Iy = Iy;  
           obj.Iz = Iz;
           obj.Ixz = Ixz;

           % Coefficients and Derivatives
           obj.CL0 = coeff.CL0;
           obj.CLalf = coeff.CLalf; 
           obj.CLq = coeff.CLq;
           obj.CLdf = coeff.CLdf;
           obj.CLde = coeff.CLde;
           obj.CD0 = coeff.CD0; 
           obj.CD1 = coeff.CD1;
           obj.CD2 = coeff.CD2;
           obj.CYb = coeff.CYb;
           obj.CYp = coeff.CYp;
           obj.CYr = coeff.CYr;
           obj.CYda = coeff.CYda;
           obj.CYdr = coeff.CYdr;
           obj.Clb = coeff.Clb;
           obj.Clp = coeff.Clp;
           obj.Clr = coeff.Clr;
           obj.Clda = coeff.Clda;
           obj.Cldr = coeff.Cldr;
           obj.Cm0 = coeff.Cm0;
           obj.Cmalf = coeff.Cmalf;
           obj.Cmq = coeff.Cmq;
           obj.Cmdf = coeff.Cmdf;
           obj.Cmde = coeff.Cmde;
           obj.Cnb = coeff.Cnb;
           obj.Cnp = coeff.Cnp;
           obj.Cnr = coeff.Cnr;
           obj.Cnda = coeff.Cnda;
           obj.Cndr = coeff.Cndr;   
           
           obj.p_exp = p_exp;
           
       end
       function [dZ,aeroCoeff] = stateDervs(obj,Z,VR)
           % Z = [u,v,w,p,q,r,phi,thet,psi,x,y,z,df,da,de,dr,CTx,CTy]
           u = Z(1); v = Z(2); w = Z(3); p = Z(4); q = Z(5); r = Z(6);
           phi = Z(7); thet = Z(8); psi = Z(9);% x = Z(10); y = Z(11); 
           z = Z(12); df = Z(13); da = Z(14); de = Z(15); dr = Z(16);
           % CTx = Z(17); CTy = Z(18); 
           CTx = 0; CTy = 0;
           dZ = zeros(12,1);
           
           %%%% Wind model
           % exponential
           Wx = VR*(-z)^obj.p_exp;
           Wxz = (obj.p_exp*VR)*((-z)^obj.p_exp)/z;
           
           %%%% 
           
           Cphi = cos(phi); Sphi = sin(phi);
           Cthet = cos(thet); Sthet = sin(thet);
           Cpsi = cos(psi); Spsi = sin(psi);
           
           dZ(10) = u*Cpsi*Cthet + v*(Cpsi*Sthet*Sphi-Spsi*Cphi) + w*(Cpsi*Sthet*Cphi+Spsi*Sphi) + Wx;
           dZ(11) = u*Spsi*Cthet + v*(Spsi*Sthet*Sphi+Cpsi*Cphi) + w*(Spsi*Sthet*Cphi-Cpsi*Sphi);
           dZ(12) = -u*Sthet + v*Cthet*Sphi + w*Cthet*Cphi;
           
           dZ(7) = p + q*Sphi*Sthet/Cthet + r*Cphi*Sthet/Cthet;
           dZ(8) = q*Cphi - r*Sphi;
           dZ(9) = q*Sphi/Cthet + r*Cphi/Cthet;
           
           VT = sqrt(u^2 + v^2 + w^2);
           aoa = atan(w/u);
           bet = asin(v/VT);
           Q = 0.5*1.225*VT*VT;
           
           Ca = cos(aoa); Sa = sin(aoa);
           Cb = cos(bet); Sb = sin(bet);
           
           % body frame ang. vel. not expressed in wind axes.
           % F1 = D ; F2 = Y ; F3 = L
           CL = ( obj.CL0 + obj.CLalf*aoa + obj.CLq*(0.5*q*obj.c/VT) + obj.CLdf*df + obj.CLde*de );
           CD = ( obj.CD0 + obj.CD1*CL + obj.CD2*CL^2 );
           CY = ( obj.CYb*bet + obj.CYp*(0.5*p*obj.b/VT) + obj.CYr*(0.5*r*obj.b/VT) + obj.CYda*da + obj.CYdr*dr );
           F3 = Q*obj.S*CL;
           F1 = Q*obj.S*CD;
           F2 = Q*obj.S*CY;
           
           Fx = -F1*Ca*Cb - F2*Ca*Sb + F3*Sa + Q*obj.S*CTx;
           Fy = -F1*Sb + F2*Cb + Q*obj.S*CTy;
           Fz = -F1*Cb*Sa - F2*Sa*Sb - F3*Ca;
           
           % M1 = L' ; M2 = M' ; M3 = N';
           Cl = ( obj.Clb*bet + obj.Clp*(0.5*p*obj.b/VT) + obj.Clr*(0.5*r*obj.b/VT) + obj.Clda*da + obj.Cldr*dr );
           Cm = ( obj.Cm0 + obj.Cmalf*aoa + obj.Cmq*(0.5*q*obj.c/VT) + obj.Cmdf*df + obj.Cmde*de );
           Cn = ( obj.Cnb*bet + obj.Cnp*(0.5*p*obj.b/VT) + obj.Cnr*(0.5*r*obj.b/VT) + obj.Cnda*da + obj.Cndr*dr );
           M1 = Q*obj.S*obj.b*Cl;
           M2 = Q*obj.S*obj.c*Cm;
           M3 = Q*obj.S*obj.b*Cn;
           
           Mx = M1*Ca*Cb - M2*Ca*Sb - M3*Sa;
           My = M1*Sb + M2*Cb;
           Mz = M1*Cb*Sa - M2*Sa*Sb + M3*Ca + Q*obj.S*obj.d*CTy;
           
           aeroCoeff = [CD,CY,CL,Cl,Cm,Cn];
           
           g = 9.806;
           dZ(1) = -w*q + v*r + Fx/obj.m - g*Sthet - Wxz*Cpsi*Cthet*dZ(12); 
           dZ(2) = -r*u + p*w + Fy/obj.m + g*Cthet*Sphi - Wxz*( Cpsi*Sthet*Sphi - Spsi*Cphi )*dZ(12);
           dZ(3) = -p*v + q*u + Fz/obj.m + g*Cthet*Cphi - Wxz*( Cpsi*Sthet*Cphi + Spsi*Sphi )*dZ(12);
           
           term1 = Mx - q*r*(obj.Iz - obj.Iy) + p*q*obj.Ixz;
           term2 = Mz - p*q*(obj.Iy - obj.Ix) - q*r*obj.Ixz;
           term3 = obj.Ix*obj.Iz - obj.Ixz^2;
           
           dZ(4) = ( obj.Iz*term1 + obj.Ixz*term2 )/term3;
           dZ(5) = ( My - p*r*(obj.Ix - obj.Iz) - (p^2 - r^2)*obj.Ixz )/obj.Iy;
           dZ(6) = ( obj.Ixz*term1 + obj.Ix*term2 )/term3;
           
       end
   end
end