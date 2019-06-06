clear all; close all;
clc;

load('eigData6DoF.mat');
%X = [u,v,w,p,q,r,Phi,Theta,Psi,x,y,z]
%U = [df,da,de,dr,CTx,CTy]
ustateind=1;
thetastateind=8;
u0=20;
pexn=zeros(N,d);
Tperiod=dsSoln.T;
magratiovalnom=0.6;
phaseratiovalnom=1.67;
id_eigE = [];
for indFE=1:length(eigE)
    if( abs(real(eigE(indFE)))<1 && (abs(imag(eigE(indFE)))<2*pi/Tperiod*0.25*N))
        for p=1:d
            pexn(:,p)=fftshift(fft(eigVec(1+(p-1)*N:p*N,indFE)));
            pexn(1,p) = pexn(1,p)/u0;
        end
        [maxval2, maxind2]=max(sum(abs(pexn),2));% adding abs values.
        dominantFE2=eigE(indFE)+sqrt(-1)*(maxind2-N/2)*2*pi/Tperiod;
        magratioval=abs((pexn(maxind2,ustateind)/u0)/pexn(maxind2,thetastateind));
        phaseratioval=angle((pexn(maxind2,ustateind)/u0)/pexn(maxind2,thetastateind));
        if(magratioval>magratiovalnom && abs(phaseratioval-phaseratiovalnom)<1*pi/180)
            id_eigE(end+1) = eigE(indFE);
            disp([eigE(indFE) magratioval phaseratioval])
        end
    end
end

figure
hold on
plot(real(eigE),imag(eigE),'xm','LineWidth',1)
plot(real(id_eigE),imag(id_eigE),'.b','MarkerSize',10)
title('Floquet Exponents');
legend('Spectral FE','Identified FE');
grid minor

%%

[grp,idx,trim_ind] = estEigGroups_temp(eigE,N,dsSoln.T);
eigVec(:,trim_ind) = [];
eigE(trim_ind) = [];

figure
hold on
plot(real(eigE),imag(eigE),'xm','LineWidth',1)
for i = 1:numel(grp)
    plot(real(grp{i}),imag(grp{i}),'.b','MarkerSize',5);
end
x = linspace(1.2*min(real(eigE)),1.2*max(real(eigE)),100);
plot(x,0.5*pi*N/Tperiod*ones(1,100),'-k');
plot(x,-0.5*pi*N/Tperiod*ones(1,100),'-k');
grid minor

pexn = cell(1,length(grp));
for i = 1:length(grp)
    for j = 1
%     for j = 1:length(idx{i})
        indFE = idx{i}(j);
        for p=1:d
            pexn{i}(:,p) = fftshift(fft( eigVec( 1+(p-1)*N:p*N , indFE ) ) );  
        end
        pexn{i}(:,1) = pexn{i}(:,1)/u0;
        [maxval2, maxind2] = max(sum(abs(pexn{i}),2));% adding abs values.
        dominantFE2 = eigE(indFE)+sqrt(-1)*(maxind2-N/2)*2*pi/Tperiod;
        magratioval=abs((pexn{i}(maxind2,ustateind))/pexn{i}(maxind2,thetastateind));
        phaseratioval=angle((pexn{i}(maxind2,ustateind))/pexn{i}(maxind2,thetastateind));
        display([i dominantFE2 magratioval phaseratioval maxind2])
%         if(magratioval>magratiovalnom && abs(phaseratioval-phaseratiovalnom)<1*pi/180)
%             id_eigE(end+1) = eigE(indFE);
%             disp([i eigE(indFE) magratioval phaseratioval])
%         end
    end
end
