function f = objectiveFunc(X, trajData, objType)
    nState = 12; nControl = 6; nOptim = nState + nControl;
    if strcmp(trajData.shape, 'loiter')
        if strcmp(trajData.type, 'not-same')
            N = (length(X)-3)/nOptim;
            VR = X(nOptim*N+2); trajData.VR = VR;
        else
            N = (length(X)-2)/nOptim;
        end
    else
        if strcmp(trajData.type, 'not-same')
            N = (length(X)-2)/nOptim;
            VR = X(nOptim*N+2); trajData.VR = VR;
        else
            N = (length(X)-1)/nOptim;
        end
    end
    
    if strcmp(objType, 'VR')
        if strcmp(trajData.type, 'same')
            error("Cannot optimize VR if wind shear kept constant!")
        else
            f = VR;
        end
    end
end