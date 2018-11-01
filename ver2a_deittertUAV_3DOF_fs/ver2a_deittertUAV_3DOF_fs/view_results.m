clear all
% load temp_initguess
load loiter_002
load acmod_001

M = (length(X)-5)/6;
N = 100;

VR = X(end);
tf = X(end-1);
t = linspace(0,tf,N);

fourCoeff(3,2*M+1) = 0;
len = 2*M+1;
for i = 1:3
    fourCoeff(i,:) = X((i-1)*len+1:i*len);
end

x(N,1) = 0; y(N,1) = 0; z(N,1) = 0;
dx(N,1) = 0; dy(N,1) = 0; dz(N,1) = 0;
ddx(N,1) = 0; ddy(N,1) = 0; ddz(N,1) = 0;

V(N,1) = 0; chi(N,1) = 0; gam(N,1) = 0;
CL(N,1) = 0; mu(N,1) = 0; CT(N,1) = 0;
for j = 1:N
    bases = fBasis(t(j),tf,M);
    x(j) = fourCoeff(1,:)*bases(1,:)';
    y(j) = fourCoeff(2,:)*bases(1,:)';
    z(j) = fourCoeff(3,:)*bases(1,:)';

    dx(j) = fourCoeff(1,:)*bases(2,:)';
    dy(j) = fourCoeff(2,:)*bases(2,:)';
    dz(j) = fourCoeff(3,:)*bases(2,:)';

    ddx(j) = fourCoeff(1,:)*bases(3,:)';
    ddy(j) = fourCoeff(2,:)*bases(3,:)';
    ddz(j) = fourCoeff(3,:)*bases(3,:)'; 

    %%% Wind model
    % exponential
    p_exp = 1;
    Wx = VR*(-z(j)).^p_exp;
    Wxz = (p_exp*VR)*((-z(j)).^p_exp)./z(j);
    wind = [Wx,Wxz];

    Z = evalZ_3DOF([x(j),y(j),z(j)],[dx(j),dy(j),dz(j);ddx(j),ddy(j),ddz(j)],wind,ac);
    V(j) = Z(1);
    chi(j) = Z(2);
    gam(j) = Z(3);
    CL(j) = Z(4);
    mu(j) = Z(5);
    CT(j) = Z(6);
end

figure
plot3(x,-y,-z,'-b');
hold on
strt = plot3(x(1),-y(1),-z(1),'*g'); % start
fin = plot3(x(end),-y(end),-z(end),'or'); % end
xlabel('x');
ylabel('y');
zlabel('z');
view(61,31)
legend([strt,fin],'start','end');
title('Trajectory');
grid on
hold off

pltmarker = '-k';
figure
subplot(2,3,1)
plot(t,V,pltmarker);
title('V');

subplot(2,3,2)
plot(t,chi,pltmarker);
title('\chi');

subplot(2,3,3)
plot(t,gam,pltmarker);
title('\gamma');

subplot(2,3,4)
plot(t,CL,pltmarker);
title('C_L');

subplot(2,3,5)
plot(t,mu,pltmarker);
title('\mu');

subplot(2,3,6)
plot(t,CT,pltmarker);
title('C_{T}');
%%
% save('xyz','x','y','z','VR','tf');