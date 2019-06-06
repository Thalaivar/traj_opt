function [FE, groupSizes, AM] = identifyFloquet(eigE, N, T)
    imRelTol = 1e-2; reRelTol = 1e-3; reAbsTol = 1e-4;
    if(0.25*N < 2*pi/T)
        error("N is not large enough, must be at least 8*pi/T")
    end
    % truncate eigenvalues to chop off spurious modes at the top
    eigE = eigE((abs(imag(eigE)) <= 0.25*N*2*pi/T));
    % sort eigenvalues according to increasing real part
    sortedEigE = sort(eigE, 'descend', 'ComparisonMethod', 'real');
    eigArrLen = length(sortedEigE);
    % identify vertical lines
    j = 1; FE = []; groupSizes = []; AM = [];% grpSizes = []; maxViolation = []; groupFE = [];
    while(j+2 <= eigArrLen)
        j1 = sortedEigE(j); j2 = sortedEigE(j+1); j3 = sortedEigE(j+2);
        % check if we are on the vertical line
        if reCond(j1, j2, j3, reRelTol, reAbsTol)
            % if on maxFE line, start checking imag parts
            currGroup = [];
            while imCond(j1,j2,j3,T,imRelTol) && reCond(j1,j2,j3,reRelTol,reAbsTol)
                currGroup(end+1) = sortedEigE(j);
                j = j + 1;
                if j + 2 > eigArrLen, break; end
                j1 = sortedEigE(j); j2 = sortedEigE(j+1); j3 = sortedEigE(j+2); 
            end
            if length(currGroup) > 0.25*N
                % if we have a non empty group check termination points
                currGroup = endCondition([j1, j2, j3], T, currGroup, imRelTol, reRelTol, reAbsTol);
                [FE,AM] = examineGroup(currGroup, FE, N, AM);
                groupSizes(end+1) = length(currGroup);
            else
                j = j + max(length(currGroup), 1);
            end
            if length(groupSizes) == 2
                break;
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
    if compareRealFE(j1, j2, reRelTol, reAbsTol)
        if compareRealFE(j2, j3, reRelTol, reAbsTol)
            isSameGroup = true;
        end
    end
end
function result = compareRealFE(x, y, reRelTol, reAbsTol)
    if(abs(real(x)) < 1e-4 && abs(real(y)) < 1e-4)
        result = abs(real(x) - real(y)) < reAbsTol;
    else
        result = (abs((real(x) - real(y))/real(y)) <= reRelTol);
    end
end
function [currGroup, indx] = endCondition(J, T, currGroup, imRelTol, reRelTol, reAbsTol)
    for i = 1:length(J)
        isSameGrp = false;
        for j = 1:length(currGroup)
            if compareRealFE(J(i), currGroup(j), reRelTol, reAbsTol)
                imagCheck = abs(mod((imag(J(i))-imag(currGroup(j)))/(2*pi/T), 1));
                imagCheck = abs(round(imagCheck) - imagCheck);
                if imagCheck < imRelTol
                    if abs(imag(currGroup(j)) > 1e-7), isSameGrp = true; end
                end
            else
                break
            end
        end
        if isSameGrp
            currGroup(end+1) = J(i);
        else
            indx = i - 1;
            break
        end
    end
end
function [FE,AM] = examineGroup(currGroup, FE, N, AM)
    groupSize = length(currGroup);
    if groupSize < 0.5*N-6 
        errmsg = ['Too few points in vertical line corresponding to : ', num2str(real(currGroup(1)))];
        error(errmsg);
    else
        % if N/2 points, it correpsonds to a real FE
        if groupSize >= 0.5*N - 6 && groupSize <= 0.5*N + 6
            % check for spectral eigenvalues on real axis
            [~, indx] = find(abs(imag(currGroup)) < 1e-7);
            AM(end+1) = 1;
            if length(indx) ~= 1
                % if no real spectral eigs, investigate
                warnmsg = ['Incorrect number of real spectral eigs for line corresponding to real FE = ',num2str(real(currGroup(1)))];
                warning(warnmsg);
                
                if length(indx) > 1
                    % in most cases, real spurious eigs close to vertical
                    % line are added to group, picking the first real
                    % eigenvalue in currGroup would usually correspond to
                    % the actual real FE
                    FE(end+1) = currGroup(indx(1));
                    currGroup = [currGroup(find(abs(imag(currGroup)) > 1e-7)), currGroup(indx(1))];
                else
                    [~, indx] = min(abs(imag(currGroup)));
                    FE(end+1) = real(currGroup(indx));
                end
            else
                % add real FE to estimate list
                FE(end+1) = currGroup(indx);
            end
        % if N points, either conjugate pair or real FE with AM = 2
        elseif groupSize >= N - 6 && groupSize <= N + 6
            [~, indx] = find(abs(imag(currGroup)) < 1e-7);
            if isempty(indx)
                % zero real spectral eigenvalues => conjugate pair FE
                [~, indx] = min(abs(imag(currGroup)));
                FE(end+1) = currGroup(indx);
                FE(end+1) = complex(real(currGroup(indx)), -imag(currGroup(indx)));
                AM(end+1) = 1;
            elseif length(indx) == 2
                % two real spectral eigenvalues => real FE with AM = 2
                FE(end+1) = currGroup(indx(1)); 
                FE(end+1) = currGroup(indx(2)); 
                AM(end+1) = 2;
            else
                warnmsg = ['Incorrect number of real spectral eigs for line corresponding to N points = ',num2str(real(currGroup(1)))];
                warning(warnmsg);
                [~,indx] = min(abs(imag(currGroup)));
                FE(end+1) = complex(real(currGroup(indx)), imag(currGroup(indx)));
                FE(end+1) = complex(real(currGroup(indx)), -imag(currGroup(indx)));
                AM(end+1) = 1;
            end
        % if 3N/2 points, either complex conjugates and one real FE, or all
        % three real with AM = 3
        elseif groupSize >= 1.5*N - 6 && groupSize <= 1.5*N + 6
            [~, indx] = find(abs(imag(currGroup)) < 1e-7);
            if length(indx) == 1
                % one real spectral eigenvalue => two complex conjugates
                % and one real FE
                FE(end+1) = currGroup(indx);
                [~, indx] = min(abs(imag(currGroup(find(abs(imag(currGroup)) > 1e-7)))));
                FE(end+1) = currGroup(indx);
                FE(end+1) = complex(real(currGroup(indx)), -imag(currGroup(indx)));
                AM(end+1) = 1; AM(end+1) = 1;
            elseif length(indx) == 3
                % three real spectral eigenvalues => real FE with AM = 3
                FE(end+1:end+3) = currGroup(indx);
                AM(end+1) = 3;
            else
                warnmsg = ['Incorrect number of real spectral eigs for line corresponding to 3N/2 points = ',num2str(real(currGroup(1)))];
                warning(warnmsg);
                [~,indx] = min(abs(imag(currGroup)));
                FE(end+1) = complex(real(currGroup(indx)), imag(currGroup));
                FE(end+1) = complex(real(currGroup(indx)), -imag(currGroup));
                FE(end+1) = real(currGroup(indx));
                AM(end+1) = 1; AM(end+1) = 1;
            end
        else
            errmsg = ['Too many points in vertical line corresponding to : ', num2str(real(currGroup(1)))];
            error(errmsg);
        end
    end
end
