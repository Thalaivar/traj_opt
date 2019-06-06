function ac = create_acmod()
    m = 4.5;
    S = 0.473;
    c = 0.1614;
    b = 3;
    d = 0.11533;
    Ix = 1.4;
    Iy = 0.9;
    Iz = 2.3;
    Ixz = 0.05;

    % Stability derivatives and coefficients (VT = 20m/s)
    coeff = struct;

    coeff.CL0 = 0.312;
    coeff.CLalf  = 5.926;
    coeff.CLq = 8.240;
    coeff.CLdf = 3.508;
    coeff.CLde = 0.384;

    coeff.CD0 = 0.01733;
    coeff.CD1 = -0.0337;
    coeff.CD2 = 0.0517;

    coeff.CYb = -0.428;
    coeff.CYp = 0.021; 
    coeff.CYr = 0.279;
    coeff.CYda = 0.131;
    coeff.CYdr = -0.391;
    coeff.Clb = -0.029;
    coeff.Clp = -0.673;
    coeff.Clr = 0.113;
    coeff.Clda = -0.691;
    coeff.Cldr = -0.021;
    coeff.Cm0 = 0.069;
    coeff.Cmalf = -0.392;
    coeff.Cmq = -36.383;
    coeff.Cmdf = 1.075;
    coeff.Cmde = -2.195;
    coeff.Cnb = 0.122;
    coeff.Cnp = -0.036;
    coeff.Cnr = -0.081;
    coeff.Cnda = -0.045;
    coeff.Cndr = 0.116;

    p_exp = 1;

    % initializing the aircraft object
    % caters for 4 CS: df da de dr
    % quadratic expansion for CD
    ac = aircraft6DoF(m,S,b,c,d,Ix,Iy,Iz,Ixz,coeff,p_exp);
end