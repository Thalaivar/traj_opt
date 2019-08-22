function [FE,PE] = getLam0PE0(eigE,eigVec,tfin,n,d,N)

n0 = 0;

pibT = pi/tfin;
t = ((tfin/N):(tfin/N):tfin).';

[grp,ind,trim_ind] = estEigGroups(eigE,N,tfin,n);
eigE(trim_ind) = [];
eigVec(:,trim_ind) = [];

n = length(grp);

% grpInd = 25; % default member index within the sorted group

omg0(n,1) = 0;
m(n,1) = 0;
FE(n,1) = 0;
PE = cell(1,n);
mInd(n,1) = 0;

tolVal = 1e-4;

for i = 1:n
    
    [~,sortInd] = sort(abs(imag(grp{i})),'descend');
    
    MM = 4;
    mInd1 = ind{i}(sortInd(end-MM));
    mInd1p1 = ind{i}(sortInd(end-MM+1));
    mInd1m1 = ind{i}(sortInd(end-MM-1));
    
    m1 = floor(imag(eigE(mInd1))/(2*pibT));
    for k = [mInd1m1,mInd1p1]
        if abs(imag(eigE(mInd1))+imag(eigE(k)))<tolVal
           m2 = floor(imag(eigE(k))/(2*pibT));
           mInd2 = k;
        end
    end
    
    omg01 = imag(eigE(mInd1))+2*pibT*(-m1);
    omg02 = imag(eigE(mInd2))+2*pibT*(-m2);
    
    if omg01<(pibT+tolVal)
        omg0(i) = omg01;
        m(i) = m1;
        mInd(i) = mInd1;
    elseif omg01>(pibT+tolVal) && omg01<(2*pibT+tolVal)
        omg0(i) = omg02;        
        m(i) = m2;
        mInd(i) = mInd2;
    else
        error('Both omg01 and omg02 do not lie within (0,2*pi/T]')
    end
    
    
    FE(i) = eigE(mInd(i)) + 1i*2*pibT*(n0-m(i));
    % FE(i) = eigE(mInd(i));
    
    PE{i} = zeros(N,d);
    for j = 1:d
        PE{i}(:,j) = eigVec((j-1)*N+1:j*N,mInd(i)).*exp( 1i*2*pibT*(-n0+m(i))*t );
        % PE{i}(:,j) = eigVec((j-1)*N+1:j*N,mInd(i));
    end
        
end
end