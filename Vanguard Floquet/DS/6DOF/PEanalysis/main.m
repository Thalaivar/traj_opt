clear all
load acmod
load eigData6DoF
addpath('./JacEval');
ac.p_exp = 1;

uavg = mean(dsSoln.u);

Xeq = zeros(18,1);
Xeq(1) = 20; Xeq(10:12) = [10;10;10];
VR = 0;
Jac = numJacEval(ac,Xeq,VR);
Jac(10:12,:) = [];
Jac(:,10:12) = [];

[eigVec,eigVal] = eig(Jac);
eigVal = diag(eigVal);

ref_ubyThet = (eigVec(1,8)/eigVec(8,8))/uavg;
ref_phs_ubyThet = angle(ref_ubyThet);
ref_mag_ubyThet = abs(ref_ubyThet);

for i=1:9
    ubyThet = (eigVec(1,i)/eigVec(8,i))/uavg;
    phs_ubyThet = angle(ubyThet);
    mag_ubyThet = abs(ubyThet);
    if (mag_ubyThet-ref_mag_ubyThet)>=0 && abs(ref_phs_ubyThet-phs_ubyThet)<=5*pi/180
        display([i, eigVal(i)])
        display([eigVec(:,i)])
    end    
end

rmpath('./JacEval');