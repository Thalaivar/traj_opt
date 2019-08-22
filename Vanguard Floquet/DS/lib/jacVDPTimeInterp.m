function A = jacVDPTimeInterp(t, trajData)
    tt = trajData.T*trajData.fourierGrid/(2*pi);
    x = interp_sinc(tt, trajData.X(:,1), t);
    y = interp_sinc(tt, trajData.X(:,2), t);
    A = jacVDP([x,y], trajData);
end