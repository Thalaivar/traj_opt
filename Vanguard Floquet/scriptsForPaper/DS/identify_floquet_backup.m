function FE = identify_floquet_backup(eigE, N, T)
    imCondTol = 1e-3;
    vLineTol = 1e-5;
    if(0.25*N < 2*pi/T)
        error("N is not large enough, must be at least 8*pi/T")
    end
    % truncate eigenvalues to chop off spurious modes at the top
    eigE = eigE(find(abs(imag(eigE)) <= 0.25*N*2*pi/T));
    % sort eigenvalues according to increasing real part
    sortedEigE = sort(eigE, 'descend', 'ComparisonMethod', 'real');
    eigGroups = {};
    % identify vertical lines
    j = 1; FE = [];
    while(j+2 <= length(sortedEigE))
        j1 = sortedEigE(j); j2 = sortedEigE(j+1); j3 = sortedEigE(j+2);
        isVerticalLine = false;
        if abs(real(j1) - real(j2)) < vLineTol
            if abs(real(j2) - real(j3)) < vLineTol
                isVerticalLine = true;
            end
        end

        isFELine = false;
        if(isVerticalLine)
            imagCheck1 = mod(imag(j1) - imag(j2), 2*pi/T);
            imagCheck2 = mod(imag(j2) - imag(j3), 2*pi/T);
            imagCheck3 = mod(imag(j3) - imag(j1), 2*pi/T);
            if(imagCheck1 < imCondTol)
                isFELine = true;
            elseif(imagCheck2 < imCondTol)
                isFELine = true;
            elseif(imagCheck3 < imCondTol)
                isFELine = true;
            end
        end

        if(isFELine && isVerticalLine)
            verticalLineSize = 1; currFE = real(j1);
            currGroups = [];
            while(abs(real(sortedEigE(j)) - currFE/verticalLineSize) < vLineTol)
                currFE = currFE + real(sortedEigE(j));
                verticalLineSize = verticalLineSize + 1;
                currGroups(end+1) = real(sortedEigE(j));
                j = j + 1;
                if j > length(sortedEigE)
                    break
                end
            end
            eigGroups{end+1} = currGroups;
            if verticalLineSize < 0.75*N && vertical
                FE(end+1) = currFE/verticalLineSize;
            else
                if verticalLineSize > 1.25*N
                    FE(end+1) = currFE/verticalLineSize;
                    FE(end+1) = currFE/verticalLineSize;
                    FE(end+1) = currFE/verticalLineSize;
                elseif verticalLineSize > 0.75*N && verticalLineSize < 1.25*N
                    FE(end+1) = currFE/verticalLineSize;
                    FE(end+1) = currFE/verticalLineSize;
                end
            end
        else
            j = j + 1;
        end
    end
end