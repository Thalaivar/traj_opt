function [mode, grpSize] = modeID(eigE, eigVec, d, N, T, floqRef, urms)
    n = 5; % no. groups to be estimated
    maxFreq = 10;
    wrappi=@(x)(x-floor(x/pi)*2*(pi));
    
    [floqTst,grp] = pickFloq(eigE,eigVec,n,d,N,T,urms,maxFreq);
    
    % identify the mode closest to a group of the reference solution
    grpIdx = 2;
    modeAnglRsd(n,1) = 0;
    modeMagRsd(n,1) = 0;
    rsdAnglMat(maxFreq,9) = 0;
    for i = 1:n
    %     rsdAnglMat = angle(floqRef.PEhat(grpIdx,:,:)./floqTst.PEhat(i,:,:));
        for k = 1:2*maxFreq
            rsdAnglMat(k,:) = wrappi(abs(angle(floqRef.PEhat(i,k,:)/floqRef.PEhat(i,k,1)) - ...
                                   angle(floqTst.PEhat(grpIdx,k,:)./floqTst.PEhat(grpIdx,k,1))));
        end
        rsdMagMat = abs(floqRef.PEhat(grpIdx,:,:)) - abs(floqTst.PEhat(i,:,:));
        modeAnglRsd(i) = norm(rsdAnglMat(:));
        modeMagRsd(i) = norm(rsdMagMat(:));
    end

    scoreAnglRsd = modeAnglRsd/norm(modeAnglRsd);
    scoreMagRsd = modeMagRsd/norm(modeMagRsd);
    scoreComb =0.5*(scoreMagRsd+ scoreAnglRsd);


    [minscoreangle, minscoreangleind]=min(scoreAnglRsd);
    [minscoremag, minscoremagind]=min(scoreMagRsd);
    [minscorecombined, minscorecombinedind]=min(scoreComb);

%     modeInd=0;%not required. Helps debug

    if((minscoreangle>.2)&&(minscoremag>0.2))
         modeInd=minscorecombinedind;
    elseif(minscoreangle<.2)
        modeInd=minscoreangleind;
    else
        modeInd=minscoremagind;   
    end

    % if(minscoreangle<.2)
    %     modeInd=minscoreangleind;
    % elseif(minscoremag<0.2)
    %     modeInd=minscoremagind;
    % else
    %     modeInd=minscorecombinedind;
    % end
    % 
    % [minScoreAngl, minScoreAnglInd] = min(scoreAnglRsd);
    % [minScoreMag, minScoreMagInd] = min(scoreMagRsd);
    % scoreComb = 0.5*(scoreMagRsd+scoreAnglRsd);
    % [minScoreComb, minScoreCombInd] = min(scoreComb);

    % modeInd=0;
    % if abs( max( scoreAnglRsd ) - sqrt(0.2) )/sqrt(0.2) < 1e-1
    % 	if(minScoreAngle < 0.1) % if phase score is very good, employ only that
    % 		modeInd = minScoreAngInd;
    % 	else
    % 		modeInd = minScoreCombInd;
    % 	end
    % else
    % 	modeInd = minScoreMagInd;
    % end
    
    mode = real(floqTst.FE(modeInd));
    grpSize = grp(modeInd);
%     mode.scoreAnglRsd = scoreAnglRsd;
%     mode.scoreMagRsd = scoreMagRsd;
%     mode.scoreComb = scoreComb;
end