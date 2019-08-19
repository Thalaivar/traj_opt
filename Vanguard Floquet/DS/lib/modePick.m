function modeInfo = modePick(modeRef,eigE,eigVec,n,d,N,T,urms,maxFreq,grpIdx)
    % last modified: 17/06/19
    % 
    % test solution quantities:
    % eigE, eigeVec, N, d, T, urms
    % 
    % maxFreq: maximum abs(frequency) up to which test and reference PEhat are
    % compared
    % 
    % determine PE corresponding to the FE with imaginary part 2*pi*n0/T in each group
    % choose omg0 in [0, pi/T]
    % return the DFT of the PE associated with the identified FE and omg0 for each FE group

    n0 = 0;

    pibT = pi/T;
    t = (T/N)*(1:N).';

    [grp,ind,trim_ind] = estEigGroups(eigE,N,T,n);
    eigE(trim_ind) = [];
    eigVec(:,trim_ind) = [];

    grpInd = 3; % default member index within the sorted group

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

        FE(i) = eigE(mInd(i)) + 1i*2*pi*(n0-m(i))/T;

        for j = 1:d
            PEhat(i,:,j) = fftshift( fft( eigVec((j-1)*N+1:j*N,mInd1).*exp(1i*2*pibT*(n0-m(i))*t) ) )/N;
        end

        % normalize uhat, vhat and what with rms of u over [0, T]
        PEhat(i,:,[1,2,3]) = PEhat(i,:,[1,2,3])/urms;

        for l = 1:N
            % normalizing states at each freq with the RMS of tbe state fourier coeff. at that freq
            nrmStates = sqrt( sum(PEhat(i,l,:).^2) );
            PEhat(i,l,:) = PEhat(i,l,:)/nrmStates;
        end    
    end

    % trim PEhat
    PEhat(:,1:N/2-maxFreq,:) = [];
    PEhat(:,end-(N/2-maxFreq-1):end,:) = [];

    modeTst = struct;
    modeTst.FE = FE;
    modeTst.PEhat = PEhat;
    modeTst.omg0 = omg0;
    modeTst.pibT = pibT;

    % identify the mode closest to a group of the reference solution
    % wrappi=@(x)(x-floor(x/pi)*2*(pi));

    % grpIdx = 2; % selected group in the reference solution
    modeAnglRsd(n,1) = 0;
    modeMagRsd(n,1) = 0;
    rsdAnglMat(maxFreq,9) = 0;
    for i = 1:n

        for k = 1:2*maxFreq
    %         rsdAnglMat(k,:) = abs( wrappi( abs( angle( modeRef.PEhat(grpIdx,k,:)/modeRef.PEhat(grpIdx,k,1) ) - ...
    %                                             angle( modeTst.PEhat(i,k,:)/modeTst.PEhat(i,k,1) ) ) ) );
            rsdAnglMat(k,:) = angle( (modeRef.PEhat(grpIdx,k,:)/modeRef.PEhat(grpIdx,k,1))./(modeTst.PEhat(i,k,:)/modeTst.PEhat(i,k,1)) );
        end

        rsdMagMat = abs(modeRef.PEhat(grpIdx,:,:)) - abs(modeTst.PEhat(i,:,:));
        modeAnglRsd(i) = norm(rsdAnglMat(:));
        modeMagRsd(i) = norm(rsdMagMat(:));
    end
    scoreAnglRsd = modeAnglRsd/norm(modeAnglRsd);
    scoreMagRsd = modeMagRsd/norm(modeMagRsd);
    scoreComb = 0.5*(scoreMagRsd + scoreAnglRsd);

    [minScoreAngl, minScoreAnglInd] = min(scoreAnglRsd);
    [minScoreMag, minScoreMagInd] = min(scoreMagRsd);
    [minScoreComb, minScoreCombInd] = min(scoreComb);

    % modeInd = []; % Not required. Helps to debug
    % if((minScoreAngl>0.2)&&(minScoreMag>0.2))
    %     modeInd = minScoreCombInd;
    % elseif(minScoreAngl<0.2)
    %     modeInd = minScoreAnglInd;
    % else
    %     modeInd = minScoreMagInd;   
    % end

    modeInd = minScoreAnglInd;

    modeInfo = struct;
    modeInfo.FE = FE(modeInd);
    modeInfo.omg0 = omg0(modeInd);
    modeInfo.modeInd = modeInd;
    modeInfo.ScoreComb = minScoreComb;
    modeInfo.ScoreMag = minScoreMag;
    modeInfo.ScoreAngl = minScoreAngl;

end