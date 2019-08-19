function [floqDat,grp] = pickFloq(eigE,eigVec,n,d,N,T,urms,maxFreq)
% last modified: 13/06/19
% determine the PE corresponding to the FE with imaginary part 2*pi*n0/T
% in each group
% return the DFT of the identified FE and the omg0 for each FE group


n0 = 1;

pibT = pi/T;
t = (T/N)*(1:N).';

[grp,ind,trim_ind] = estEigGroups(eigE,N,T,n);
eigE(trim_ind) = [];
eigVec(:,trim_ind) = [];

grpInd = 10; % default member index within the sorted group

omg0(n,1) = 0;
m(n,1) = 0;
FE(n,1) = 0;
PEhat(n,N,d) = 0;
mInd(n,1) = 0;

tolVal = 1e-4;

for i = 1:n  
    
    [~,sortInd] = sort(abs(imag(grp{i})),'descend');
    
    mInd1 = ind{i}(sortInd(grpInd));
    mInd1p1 = ind{i}(sortInd(grpInd+1));
    mInd1m1 = ind{i}(sortInd(grpInd-1));
    
    m1 = floor(imag(eigE(mInd1))/(2*pibT));
    for k = [mInd1m1,mInd1p1]
        if abs(imag(eigE(mInd1))+imag(eigE(k)))<tolVal
           m2 = floor(imag(eigE(k))/(2*pibT));
           mInd2 = k;
        end
    end
    
%     omg01 = imag(eigE(mInd1))+2*pibT*(n0-m1);
%     omg02 = imag(eigE(mInd2))+2*pibT*(n0-m2);
    
    omg01 = imag(eigE(mInd1))+2*pibT*(-m1);
    omg02 = imag(eigE(mInd2))+2*pibT*(-m2);
    
    if omg01<pibT+tolVal
        omg0(i) = omg01;
        m(i) = m1;
        mInd(i) = mInd1;
    else
        omg0(i) = omg02;        
        m(i) = m2;
        mInd(i) = mInd2;
    end
    
    FE(i) = eigE(mInd(i)) + 1i*2*pi*(n0-m(i))/T;
    
    for j = 1:d
        PEhat(i,:,j) = fftshift( fft( eigVec((j-1)*N+1:j*N,mInd1).*exp(1i*2*pibT*(n0-m(i))*t) ) )/N;
    end
    
    PEhat(i,:,[1,2,3]) = PEhat(i,:,[1,2,3])/urms;

    for l = 1:N
        nrmStates = sqrt( sum(PEhat(i,l,:).^2) );
        PEhat(i,l,:) = PEhat(i,l,:)/nrmStates;
    end    
        
end

% trim PEhat
PEhat(:,1:N/2-maxFreq,:) = [];
PEhat(:,end-(N/2-maxFreq-1):end,:) = [];

floqDat = struct;
floqDat.FE = FE;
floqDat.PEhat = PEhat;
floqDat.omg0 = omg0;
floqDat.pibT = pibT;


end