clear all
load temp_initguess
load acmod_001

N = (length(X)-2)/3;
VR = X(end);
tf = X(end-1);
t = linspace(tf/N,tf,N); % note this, t = 0 not included. Periodicity is implicit

% differentiation matrices
column1 = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*pi/N)]';
D = toeplitz(column1,column1([1, N:-1:2])); % first derivative matrix
column2 = [-(N^2)/12-1/6, -((-1).^(1:(N-1)))./(2*(sin((1:(N-1))*pi/N)).^2)];
DD = toeplitz(column2,column2([1, N:-1:2])); % second derivative matrix

fac = 2*pi/tf;
    
x = X((1-1)*N+1:1*N); dx = fac*D*x; ddx = fac*fac*DD*x; 
y = X((2-1)*N+1:2*N); dy = fac*D*y; ddy = fac*fac*DD*y;   
z = X((3-1)*N+1:3*N); dz = fac*D*z; ddz = fac*fac*DD*z;

V(N,1) = 0; chi(N,1) = 0; gam(N,1) = 0;
CL(N,1) = 0; mu(N,1) = 0; CT(N,1) = 0;
for j = 1:N
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
