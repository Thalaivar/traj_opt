function [grp,idx,trim_ind] = estEigGroups_temp(eigE,N,T)
% identify as many floquet exponent groups as possible in descending order
% of the magnitude of real parts

global chkTol chkRelTol nearZeroTol
chkTol = 1e-2; % imModulo
chkRelTol = 1e-3;
nearZeroTol = 1e-3;

% trimming FEs with magnitude of imaginary parts greater than (N/4)*(2*pi/T)
    fac = 2*pi/T;
    trim_ind = (abs(imag(eigE)) > 1.01*fac*(0.25*N));
    eigE( trim_ind ) = []; % leeway of 1% is provided
    [eigE,eigInd] = sort(eigE,'descend', 'ComparisonMethod', 'real'); % sorting the eigenvalues

    
    lagVal = 5;
    
    NN = numel(eigE);

    grp = cell(1); idx = cell(1);
    grp{1} = [];  idx{1} = [];
    isGroupNext = true; % flag that determines whether another FE group exists

    k = 1; 
    grpCount = 1;
    while(isGroupNext)
        
        isGroupNext = false;
        
        while(k+lagVal-1 <= NN)
            
            if imModulo(eigE(k:k+lagVal-1),fac) && reCheck(eigE(k:k+lagVal-1))
                if isempty(grp{grpCount})
                    grp{grpCount}(end+1) = eigE(k);
                    idx{grpCount}(end+1) = eigInd(k);
                    currInd = k;
                else
                    if currGroup(grp{grpCount}(end),eigE(k))
                        grp{grpCount}(end+1) = eigE(k);
                        idx{grpCount}(end+1) = eigInd(k);
                        currInd = k;
                    else
                        isGroupNext = true;
                        grp{grpCount}(end+1:end+lagVal-1) = eigE(currInd+1:currInd+lagVal-1); % when another group exists
                        idx{grpCount}(end+1:end+lagVal-1) = eigInd(currInd+1:currInd+lagVal-1);
                        break
                    end    
                end                 
            end
            k = k+1;
        end
        if ~isGroupNext
           grp{grpCount}(end+1:end+lagVal-1) = eigE(currInd+1:currInd+lagVal-1); % when another group doesn't exists 
           idx{grpCount}(end+1:end+lagVal-1) = eigInd(currInd+1:currInd+lagVal-1);
        else
            grpCount = grpCount+1;
            grp{end+1} = [];
            idx{end+1} = [];
        end
        
    end    

end

function flagVal = currGroup(x,y)

    global chkRelTol nearZeroTol

    flagVal = false;
        if abs(real(x))<nearZeroTol
            if abs(real(y))<nearZeroTol
                flagVal = true;
            end
        else
            if abs( real(x-y)/real(x) ) < chkRelTol
                flagVal = true;
            end        
        end
        
end

function flagVal = reCheck(x)
% test if the elements of x-x(1) lie within the interval [-xchk, xchk]

    global chkRelTol nearZeroTol
    
    N = numel(x);
    
    xchk = false(1,N-1);
    
    if abs(real(x(1)))<nearZeroTol 
        for i = 1:N-1
            if abs(real( x(i+1)) ) < nearZeroTol
                xchk(i) = true;
            end
        end
    else
        for i = 1:N-1
            if abs( real( x(1)-x(i+1) )/real( x(1) ) )  < chkRelTol
               xchk(i) = true; 
            end    
        end        
    end    
    flagVal = all(xchk);
    
end

function flagVal = imModulo(x,y)

% test if the differences of elements of x-x(1) are multiples of y
    global chkTol % tolerance for mod=0

    N = numel(x);
    if N<3
        error('N>=3 is required');
    end

    counter = 0;
    for i = 2:N
        z = mod( imag(abs(x(1)) - abs(x(i))) ,y);
        if abs(z) < chkTol || abs(z-y)<chkTol
           counter = counter + 1;
        end
    end

    if counter>=2
       flagVal = true;
    else
       flagVal = false;
    end   
    
end