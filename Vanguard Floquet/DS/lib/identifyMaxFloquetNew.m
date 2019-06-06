function [maxFE, groupSize, zeroImagFE] = identifyMaxFloquetNew(eigE, N, T)
    imRelTol = 1e-4; reRelTol = 1e-2; reAbsTol = 1e-5;
    if(0.25*N < 2*pi/T)
        error("N is not large enough, must be at least 8*pi/T")
    end
    % truncate eigenvalues to chop off spurious modes at the top
    eigE = eigE((abs(imag(eigE)) <= 0.25*N*2*pi/T));
    % sort eigenvalues according to increasing real part
    sortedEigE = sort(eigE, 'descend', 'ComparisonMethod', 'real');
    eigArrLen = length(sortedEigE);
    j = 1;
    while(j+2 <= eigArrLen)
        j1 = sortedEigE(j); j2 = sortedEigE(j+1); j3 = sortedEigE(j+2);
        % check if we are on the vertical line
        if reCond(j1, j2, j3, reRelTol, reAbsTol)
            % if on maxFE line, start checking imag parts
            currGroup = [];
            while imCond(j1,j2,j3,T, imRelTol) && reCond(j1,j2,j3, reRelTol, reAbsTol)
                currGroup(end+1) = sortedEigE(j);
                j = j + 1;
                j1 = sortedEigE(j); j2 = sortedEigE(j+1); j3 = sortedEigE(j+2);
            end
            % once we are out of the loop we must check which j_i (j1, j2, j3 still ambiguous) belong to
            % the max FE line
            if ~isempty(currGroup)
                currGroup = checkJi([j1, j2, j3], T, currGroup, imRelTol, reRelTol, reAbsTol);
                groupSize = length(currGroup);
                [~,indx] = find(abs(imag(currGroup)) <= 1e-7);
                zeroImagFE = length(indx);

                if zeroImagFE > 0
                    maxFE = currGroup(indx(1));
                else
                    % if there are no real spectral eigs, we have conj pair
                    [~, indx] = min(abs(imag(currGroup)));
                    maxFE = real(currGroup(indx));
                    zeroImagFE = 1;
                end
                break
            else
                j = j + 1;
            end
        else
            j = j + 1;
        end
    end
end

function isFELine = imCond(j1, j2, j3, T, imRelTol)
    isFELine = false;
    imagCheck1 = abs(mod((imag(j1)-imag(j2))/(2*pi/T), 1));
    imagCheck2 = abs(mod((imag(j2)-imag(j3))/(2*pi/T), 1));
    imagCheck3 = abs(mod((imag(j3)-imag(j1))/(2*pi/T), 1));
    
    imagCheck1 = abs(round(imagCheck1) - imagCheck1);
    imagCheck2 = abs(round(imagCheck2) - imagCheck2);
    imagCheck3 = abs(round(imagCheck3) - imagCheck3);
    
    if(imagCheck1 < imRelTol)
        isFELine = true;
    elseif(imagCheck2 < imRelTol)
        isFELine = true;
    elseif(imagCheck3 < imRelTol)
        isFELine = true;
    end
end
function isSameGroup = reCond(j1, j2, j3, reRelTol, reAbsTol)
    isSameGroup = false;
    if compareFE(j1, j2, reRelTol, reAbsTol)
        if compareFE(j2, j3, reRelTol, reAbsTol)
            isSameGroup = true;
        end
    end
end
function result = compareFE(x, y, reRelTol, reAbsTol)
    if(abs(real(x)) < 1e-5 && abs(real(y)) < 1e-5)
        result = abs(real(x) - real(y)) <= reAbsTol;
    else
        result = (round(abs((real(x) - real(y))/real(y)),6) <= reRelTol);
    end
end
function currGroup = checkJi(Ji, T, currGroup, imRelTol, reRelTol, reAbsTol)
    for i = 1:length(Ji)
        isSameGroup = false;
        for j = 1:length(currGroup)-1
            if compareFE(currGroup(j), Ji(i), reRelTol, reAbsTol)
                 imagCheck = abs(mod((imag(Ji(i))-imag(currGroup(j)))/(2*pi/T), 1));
                 imagCheck = abs(round(imagCheck) - imagCheck);
                 if imagCheck < imRelTol
                     if abs(imag(currGroup(j))) > 1e-7
                        isSameGroup = true;
                     end
                 end
            else
                break
            end
        end
        if isSameGroup
            currGroup(end+1) = Ji(i);
        end
    end
end