function FE = identify_floquet(eigE, N, T)
    tol = 1e-4;
    if(0.25*N < 2*pi/T)
        error("N is not large enough, must be at least 8*pi/T")
    end
    % truncate eigenvalues to chop off spurious modes at the top
    eigE = eigE(find(abs(imag(eigE)) <= (0.25*N + 1)*2*pi/T));
    % sort eigenvalues according to increasing real part
    sorted_eigE = sort(eigE, 'ComparisonMethod', 'real');
    % identify vertical lines
    j = 1;
    FE = [];    
    while(j+2 <= length(sorted_eigE))
        j1 = sorted_eigE(j); j2 = sorted_eigE(j+1); j3 = sorted_eigE(j+2);
%         [real(j1), real(j2), real(j3)]
        % check if we are in a group of points belonging to vertical line
        vertcheck = false;
        if(abs(real(j1) - real(j2)) <= tol)
            if(abs(real(j2) - real(j3)) <= tol)
                vertcheck = true;
            end
        end
        % if we are in vertical line, perform check for variation in imag
        if(vertcheck)
            imagcheck = false;
            imag_diff1 = abs(imag(j1) - imag(j2))/(2*pi/T);
            imag_diff2 = abs(imag(j2) - imag(j3))/(2*pi/T); 
            imag_diff3 = abs(imag(j3) - imag(j1))/(2*pi/T); 
            if(abs(imag_diff1 - round(imag_diff1)) <= tol)
                imagcheck = true;
            elseif(abs(imag_diff2 - round(imag_diff2)) <= tol)
                imagcheck = true;
            elseif(abs(imag_diff3 - round(imag_diff3)) <= tol)
                imagcheck = true;
            end
            if(imagcheck)
                FE(end+1) = real(j1);
            end
        end
        % go to next group of eigE
        while(abs(real(sorted_eigE(j)) - real(j1)) <= tol)
            j = j + 1;
            if(j+3 > length(sorted_eigE))
                break
            end
        end
    end
end