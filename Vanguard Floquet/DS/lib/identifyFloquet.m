function FE = identifyFloquet(eigE, N, T)
    reRelTol = 1e-3;
    if(0.25*N < 2*pi/T)
        error("N is not large enough, must be at least 8*pi/T")
    end
    % truncate eigenvalues to chop off spurious modes at the top
    eigE = eigE((abs(imag(eigE)) <= 0.25*N*2*pi/T));
    % sort eigenvalues according to increasing real part
    sortedEigE = sort(eigE, 'descend', 'ComparisonMethod', 'real');
    eigArrLen = length(sortedEigE);
    % identify vertical lines
    j = 1; FE = []; grpSizes = []; maxViolation = [];
    while(j+2 <= eigArrLen)
        j1 = sortedEigE(j); j2 = sortedEigE(j+1); j3 = sortedEigE(j+2);
        isFELine = imCond(j1, j2, j3, T);
        
        if(isFELine)
            isSameGroup = reCond(j1, j2, j3, reRelTol);
            if(isSameGroup)
                groupSize = 0; meanFE = 0; currGroup = [];
                while(abs((real(sortedEigE(j)) - real(j1))) < 1e-4)
                    groupSize = groupSize + 1;
                    meanFE = meanFE + real(sortedEigE(j));
                    currGroup(end+1) = sortedEigE(j);
                    j = j + 1;
                    if j + 1 >= length(sortedEigE)
                        break
                    end
                end
                grpSizes(end+1) = groupSize;
                FE(end+1) = meanFE/groupSize;
                maxViolation(end+1) = abs(real((currGroup(1)-currGroup(groupSize)))/real(j1));
            else
                j = j + 1;
            end
        else
            j = j + 1;
        end
    end
%     FE = [FE; grpSizes; maxViolation];
end

function isFELine = imCond(j1, j2, j3, T)
    imRelTol = 1e-4;
    isFELine = false;
    imagCheck1 = abs(1 - mod((imag(j1)-imag(j2))/(2*pi/T), 1));
    imagCheck2 = abs(1 - mod((imag(j2)-imag(j3))/(2*pi/T), 1));
    imagCheck3 = abs(1 - mod((imag(j3)-imag(j1))/(2*pi/T), 1));
    if(imagCheck1 < imRelTol)
        isFELine = true;
    elseif(imagCheck2 < imRelTol)
        isFELine = true;
    elseif(imagCheck3 < imRelTol)
        isFELine = true;
    end
end
function isSameGroup = reCond(j1, j2, j3, reRelTol)
    isSameGroup = false;
    if abs((real(j1) - real(j2))/real(j1)) < reRelTol
        if abs((real(j2) - real(j3))/real(j1)) < reRelTol
            isSameGroup = true;
        end
    end
end
function FE = checkConjugate(