clear all
load acmod_001;
load loiter_001_fourCoeff
tf = 11.9998;
VR = 0.0931;
M = 8;
X(3*(2*M+1)+2,1) = 0;
X(end-1) = tf;
X(end) = VR;
for i = 1:3
    X((i-1)*(2*M+1)+1:i*(2*M+1)) = fourCoeff(i,:);
end
N = 100;
[c,ceq,dc,dceq] = Cfun(X,ac,M,N);
figure
% subplot(1,2,1)
plot(c)
title('inequality constraints');
% subplot(1,2,2)
% plot(ceq);
