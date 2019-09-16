function [FE, eigE, AM, groupSizes, eigVec] = spectralMethod(trajData)
    % change this parameter according to model
    d = 4;
    % spectral method runs on same number of grid points as trajectory
    D = trajData.D; N = trajData.N; T = trajData.T;
    
    Dmat = zeros(d*N); Mmat = zeros(d*N);
    for i = 1:N
        A = sysModel(trajData.p, trajData.X(i,:), [trajData.U(i,:), trajData.VR]);
%         A = sysModel(trajData.ac, trajData.X(i,:), trajData.U(i,:), trajData.VR);
        for j = 1:d
           for k = 1:d
                Mmat(i+(j-1)*N,(k-1)*N+i) = A(j,k); 
           end
           Dmat((j-1)*N+1:j*N,(j-1)*N+1:j*N) = D*(2*pi/T);
        end
    end
    
    [eigVec, eigE] = eig(Dmat - Mmat);
    eigE = -1*diag(eigE);
    
    [FE,groupSizes,AM] = identifyFloquet(eigE, N, trajData.T);
%     [grp,idx,trim_ind] = estEigGroups(eigE,N,T,[]);
end

function A = sysModel(p, Z, U)
    prm.m = 4.5;
    prm.S = 0.473;
    prm.CD0 = 0.0173;
    prm.CD1 = -0.0337;
    prm.CD2 = 0.0517;
    prm.p_exp = p;
    % evaluate jacobian
    A = JacEval(Z, U, prm);
    A = [A(1:3,1:3), A(1:3,6); A(6,1:3), A(6,6)];
%     A = A(1:3,1:3);
end

% function A = sysModel(ac, Z, U, VR)
%     % Z = [u,v,w,p,q,r,phi,theta,psi,x,y,z]
%     % U = [df, da, de, dr]
%     u = Z(1); v = Z(2); w = Z(3); p = Z(4); q = Z(5); r = Z(6); 
%     phi = Z(7); theta = Z(8); psi = Z(9); x = Z(10); y = Z(11); z = Z(12);
%     df = U(1); da = U(2); de = U(3); dr = U(4); CTx = 0; CTy = 0;
% 
%     Xval = [u,v,w,p,q,r,phi,theta,psi,x,y,z,df,da,de,dr,CTx,CTy].';
% 
%     A = numJacEval(ac,Xval,VR); % finite-diff Jacobian
% 
%     % d=9
%     A(:,10:11) = [];
%     A(10:11,:) = [];
% end

% function A = vanderpol_jac(t, tt, X)
%     mu = 0.4;
%     N = 0.5*(length(X)-1);
%     x = interp_sinc(tt, X(1:N,1), t);
%     xdot = interp_sinc(tt, X(N+1:2*N,1), t);
%     A = zeros(2);
%     A(1,2) = 1;
%     A(2,1) = -1 - 2*mu*x*xdot; A(2,2) = mu*(1-x^2);
% end