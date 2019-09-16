function A = jacGlycolyticTimeInterp(t, trajData)
    tt = trajData.T*trajData.fourierGrid/(2*pi);
    x = interp_sinc(tt, trajData.X(:,1), t);
    y = interp_sinc(tt, trajData.X(:,2), t);
    A = jacGlycolytic([x,y], trajData);
end