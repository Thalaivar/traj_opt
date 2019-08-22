function FE = getLam0(grp,tfin)

n0 = 0;

pibT = pi/tfin;

n = length(grp);
grpInd = 3; % default member index within the sorted group
omg0(n,1) = 0;
m(n,1) = 0;
FE(n,1) = 0;
mInd(n,1) = 0;

tolVal = 1e-4;

for i = 1:n
    [~,sortInd] = sort(abs(imag(grp{i})),'descend');

    mInd1 = sortInd(grpInd);
    mInd1p1 = sortInd(grpInd+1);
    mInd1m1 = sortInd(grpInd-1);

    m1 = floor(imag(grp{i}(mInd1))/(2*pibT));
    for k = [mInd1m1,mInd1p1]
        if abs(imag(grp{i}(mInd1))+imag(grp{i}(k)))<tolVal
           m2 = floor(imag(grp{i}(k))/(2*pibT));
           mInd2 = k;
        end
    end

    omg01 = imag(grp{i}(mInd1))+2*pibT*(-m1);
    omg02 = imag(grp{i}(mInd2))+2*pibT*(-m2);

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


    FE(i) = grp{i}(mInd(i)) + 1i*2*pibT*(n0-m(i));
end
end