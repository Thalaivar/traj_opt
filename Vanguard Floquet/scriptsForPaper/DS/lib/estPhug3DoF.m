function phugMode = estPhug3DoF(eigE,eigVec,ind,sRef)
% last modified: 9:52, 08/06/19
% for 3DoF aircraft
% eigE and eigVec timmed using the output from estEigGroups
% ind maps the elements of a group with the corresponding location in trimmed eigE
phugMode = [];

n = length(ind);
d = sRef.d;
N = sRef.N;
T = sRef.T;
Vrms = sRef.Vrms;
phugRef = sRef.phugRef;

for i = 1:n
    grpInd = ind{i}(1);
    
    % checking if the element picked from the group is real
    j = 1;
    while abs( imag( eigE(grpInd) ) ) < 1e-6
        j = j+1;
        grpInd = ind{i}(j);
    end
    
    
    PEhat = zeros(N,d);
    for k = 1:d
        PEhat(:,k) = fftshift( fft( eigVec((k-1)*N+1:k*N,grpInd) ) );
    end
    PEhat(:,1) = PEhat(:,1)/Vrms;
    
    [~, domFreq] = max( sum( abs( PEhat(:,[1,3]) ) , 2 ) );
    domFE = eigE(grpInd) + sqrt(-1)*( domFreq-N/2-1 )*2*pi/T;
    
    phugPar = PEhat(domFreq,1)/PEhat(domFreq,3);
    phugParMag = abs(phugPar);
    phugParPhs = angle(phugPar);
    
    if abs(real(domFE))>1e-4 % skip if Re(FE)~0
        if abs(PEhat(domFreq,1))>1e-4 % skip if |Vhat/V0|~0 at dominant frequency
            if ~(abs(mod(imag(domFE),2*pi/T))<1e-4 || abs(mod(imag(domFE),2*pi/T)-2*pi/T)<1e-4) % skip if the omg_0=0. Assume that Phugoid mode for the current geometry will have omg_0~=0
                if (phugParMag > phugRef.mag) && (abs( abs(phugParPhs) - abs(phugRef.phs) ) < 10*pi/180) % Phugoid condition
                    phugMode = struct;
                    phugMode.domFE = domFE;
                    phugMode.parMag = phugParMag;
                    phugMode.parPhs = phugParPhs;
                end
            end
        end
    end
    
end

    if isempty(phugMode)
        "CHeck!"
    end
end