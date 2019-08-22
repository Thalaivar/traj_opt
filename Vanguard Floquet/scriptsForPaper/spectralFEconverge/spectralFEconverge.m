load rawMaterial\acmod.mat
ac.p_exp = 1;

load('rawMaterial\result_N400_OLin_6DoF_purna.mat');

clear('eigE','eigVec','N','tfin','VR','urms')

NN = linspace(40, 400, 37);
% NN = 400;
FE_LMD(3,length(NN)) = 0;
tfin = Z(end);
N = (length(Z)-3)/16;
grps = cell(1,length(NN));

for l = 1:length(NN)
    
    D = fourierdiff(NN(l));
    t = (tfin/NN(l)):(tfin/NN(l)):tfin;    
    Dmat = zeros(d*NN(l));
    Mmat = zeros(d*NN(l));
    for i = 1:NN(l)
        A = sysModel(t(i),Z,N,ac);
        for j = 1:d
           for k = 1:d
                Mmat(i+(j-1)*NN(l),(k-1)*NN(l)+i) = A(j,k); 
           end
           Dmat((j-1)*NN(l)+1:j*NN(l),(j-1)*NN(l)+1:j*NN(l)) = D*(2*pi/tfin);
        end
    end
    [eigVec,eigE] = eig(Mmat-Dmat);  
    grp = estEigGroups(diag(eigE),NN(l),tfin,[]);
    grps{l} = grp;
    
end

save('rawMaterial\FEdat.mat','grps','NN','tfin')

figure
hold on
plot(diag(eigE),'xm')
grid minor
grpsInd = length(NN);
for i = 1:length(grps{grpsInd})
    plot(grps{grpsInd}{i},'.b')
end
title(horzcat('N = ',num2str(NN(grpsInd))));

function A = sysModel(t,Z,N,ac)

tfin = Z(end);
tgrid = (tfin/N):(tfin/N):tfin;

u = interp_sinc(tgrid,Z(0*N+1:1*N),t);
v = interp_sinc(tgrid,Z(1*N+1:2*N),t);
w = interp_sinc(tgrid,Z(2*N+1:3*N),t);
p = interp_sinc(tgrid,Z(3*N+1:4*N),t);
q = interp_sinc(tgrid,Z(4*N+1:5*N),t);
r = interp_sinc(tgrid,Z(5*N+1:6*N),t);
Phi = interp_sinc(tgrid,Z(6*N+1:7*N),t);
Thet = interp_sinc(tgrid,Z(7*N+1:8*N),t);
Psi = interp_sinc(tgrid,Z(8*N+1:9*N),t) + Z(end-2)*t;
x = interp_sinc(tgrid,Z(9*N+1:10*N),t);
y = interp_sinc(tgrid,Z(10*N+1:11*N),t);
z = interp_sinc(tgrid,Z(11*N+1:12*N),t);
df = interp_sinc(tgrid,Z(12*N+1:13*N),t);
da = interp_sinc(tgrid,Z(13*N+1:14*N),t);
de = interp_sinc(tgrid,Z(14*N+1:15*N),t);
dr = interp_sinc(tgrid,Z(15*N+1:16*N),t);
CTx = 0;
CTy = 0;


Xval = [u,v,w,p,q,r,Phi,Thet,Psi,x,y,z,df,da,de,dr,CTx,CTy].';

A = numJacEval(ac,Xval,Z(end-1)); % finite-diff Jacobian

% d=10
% A(:,10:11) = [];
% A(10:11,:) = [];

% d=9
A(:,10:12) = [];
A(10:12,:) = [];

end