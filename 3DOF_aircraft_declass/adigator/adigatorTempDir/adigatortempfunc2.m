function [adigatorFunInfo, adigatorOutputs] = adigatortempfunc2(adigatorFunInfo,adigatorInputs)
[flag, adigatorFunInfo, adigatorInputs] = adigatorFunctionInitialize(2,adigatorFunInfo,adigatorInputs);
if flag; adigatorOutputs = adigatorInputs; return; end;
solution = adigatorInputs{1};
t = adigatorInputs{2};
type = adigatorInputs{3};
nargin = 3; nargout = 1;adigatorVarAnalyzer('% model constants');
m = 4.5;
m = adigatorVarAnalyzer('m = 4.5;',m,'m',0);
rho = 1.225;
rho = adigatorVarAnalyzer('rho = 1.225;',rho,'rho',0);
S = 0.473;
S = adigatorVarAnalyzer('S = 0.473;',S,'S',0);
g = 9.806;
g = adigatorVarAnalyzer('g = 9.806;',g,'g',0);
Cd0 = 0.0173;
Cd0 = adigatorVarAnalyzer('Cd0 = 0.0173;',Cd0,'Cd0',0);
Cd1 = -0.0337;
Cd1 = adigatorVarAnalyzer('Cd1 = -0.0337;',Cd1,'Cd1',0);
Cd2 = 0.0517;
Cd2 = adigatorVarAnalyzer('Cd2 = 0.0517;',Cd2,'Cd2',0);
tf = solution.tf;
tf = adigatorVarAnalyzer('tf = solution.tf;',tf,'tf',0);
VR = solution.VR;
VR = adigatorVarAnalyzer('VR = solution.VR;',VR,'VR',0);
coeffs = solution.coeffs;
coeffs = adigatorVarAnalyzer('coeffs = solution.coeffs;',coeffs,'coeffs',0);
N = solution.N;
N = adigatorVarAnalyzer('N = solution.N;',N,'N',0);
% ADiGator IF Statement #1: START
cadaconditional1 = strcmp(type, 'analytic');
cadaconditional1 = adigatorVarAnalyzer('cadaconditional1 = strcmp(type, ''analytic'');',cadaconditional1,'cadaconditional1',0);
cadaconditional2 = strcmp(type, 'FD');
cadaconditional2 = adigatorVarAnalyzer('cadaconditional2 = strcmp(type, ''FD'');',cadaconditional2,'cadaconditional2',0);
adigatorIfInitialize(1,cadaconditional1,cadaconditional2);
adigatorIfIterStart(1,1);
    % Call to User Function get_traj --- (FunID 5)
    cadainput5_1 = t;
    cadainput5_1 = adigatorVarAnalyzer('cadainput5_1 = t;',cadainput5_1,'cadainput5_1',0);
    cadainput5_2 = tf;
    cadainput5_2 = adigatorVarAnalyzer('cadainput5_2 = tf;',cadainput5_2,'cadainput5_2',0);
    cadainput5_3 = coeffs;
    cadainput5_3 = adigatorVarAnalyzer('cadainput5_3 = coeffs;',cadainput5_3,'cadainput5_3',0);
    cadainput5_4 = N;
    cadainput5_4 = adigatorVarAnalyzer('cadainput5_4 = N;',cadainput5_4,'cadainput5_4',0);
    adigatorInputs = {cadainput5_1;cadainput5_2;cadainput5_3;cadainput5_4};
    [adigatorFunInfo, adigatorOutputs] = adigatortempfunc5(adigatorFunInfo,adigatorInputs);
    cadaoutput5_1 = adigatorOutputs{1};
    sigma = cadaoutput5_1;
    sigma = adigatorVarAnalyzer('sigma = cadaoutput5_1;',sigma,'sigma',0);
    % Call to User Function get_xu --- (FunID 3)
    cadainput3_1 = sigma;
    cadainput3_1 = adigatorVarAnalyzer('cadainput3_1 = sigma;',cadainput3_1,'cadainput3_1',0);
    cadainput3_2 = VR;
    cadainput3_2 = adigatorVarAnalyzer('cadainput3_2 = VR;',cadainput3_2,'cadainput3_2',0);
    adigatorInputs = {cadainput3_1;cadainput3_2};
    [adigatorFunInfo, adigatorOutputs] = adigatortempfunc3(adigatorFunInfo,adigatorInputs);
    cadaoutput3_1 = adigatorOutputs{1};
    X = cadaoutput3_1;
    X = adigatorVarAnalyzer('X = cadaoutput3_1;',X,'X',0);
    cadaoutput3_2 = adigatorOutputs{2};
    u = cadaoutput3_2;
    u = adigatorVarAnalyzer('u = cadaoutput3_2;',u,'u',0);
    V = X(1);
    V = adigatorVarAnalyzer('V = X(1);',V,'V',0);
    gamma = X(2);
    gamma = adigatorVarAnalyzer('gamma = X(2);',gamma,'gamma',0);
    chi = X(3);
    chi = adigatorVarAnalyzer('chi = X(3);',chi,'chi',0);
    Cl = u(1);
    Cl = adigatorVarAnalyzer('Cl = u(1);',Cl,'Cl',0);
    nu = u(2);
    nu = adigatorVarAnalyzer('nu = u(2);',nu,'nu',0);
    CT = u(3);
    CT = adigatorVarAnalyzer('CT = u(3);',CT,'CT',0);
    adigatorVarAnalyzer('% wind model');
    p_exp = 1;
    p_exp = adigatorVarAnalyzer('p_exp = 1;',p_exp,'p_exp',0);
    z = sigma(3);
    z = adigatorVarAnalyzer('z = sigma(3);',z,'z',0);
    zdot = sigma(6);
    zdot = adigatorVarAnalyzer('zdot = sigma(6);',zdot,'zdot',0);
    adigatorVarAnalyzer('%Wx = VR*(-z)^p_exp;');
    Wxz = -1*(p_exp*VR)*((-z)^(p_exp-1));
    Wxz = adigatorVarAnalyzer('Wxz = -1*(p_exp*VR)*((-z)^(p_exp-1));',Wxz,'Wxz',0);
    adigatorVarAnalyzer('% forces and drag polar');
    Cd = Cd0 + Cd1*Cl + Cd2*Cl^2;
    Cd = adigatorVarAnalyzer('Cd = Cd0 + Cd1*Cl + Cd2*Cl^2;',Cd,'Cd',0);
    L = 0.5*S*rho*Cl*V^2;
    L = adigatorVarAnalyzer('L = 0.5*S*rho*Cl*V^2;',L,'L',0);
    adigatorVarAnalyzer('% constructing the Jacobian');
    f11 = (CT - Cd)*rho*S*V/m;
    f11 = adigatorVarAnalyzer('f11 = (CT - Cd)*rho*S*V/m;',f11,'f11',0);
    f12 = Wxz*zdot*sin(chi)*cos(gamma);
    f12 = adigatorVarAnalyzer('f12 = Wxz*zdot*sin(chi)*cos(gamma);',f12,'f12',0);
    f13 = -1*g*cos(gamma) + Wxz*zdot*cos(chi)*sin(gamma);
    f13 = adigatorVarAnalyzer('f13 = -1*g*cos(gamma) + Wxz*zdot*cos(chi)*sin(gamma);',f13,'f13',0);
    f14 = 0;
    f14 = adigatorVarAnalyzer('f14 = 0;',f14,'f14',0);
    f15 = 0;
    f15 = adigatorVarAnalyzer('f15 = 0;',f15,'f15',0);
    f16 = p_exp*(p_exp-1)*VR*((-z)^(p_exp-2))*zdot*cos(chi)*cos(gamma);
    f16 = adigatorVarAnalyzer('f16 = p_exp*(p_exp-1)*VR*((-z)^(p_exp-2))*zdot*cos(chi)*cos(gamma);',f16,'f16',0);
    f21 = 0.5*rho*S*Cl*sin(nu)/(m*cos(gamma)) - (Wxz*zdot*sin(chi)/(cos(gamma)*V^2));
    f21 = adigatorVarAnalyzer('f21 = 0.5*rho*S*Cl*sin(nu)/(m*cos(gamma)) - (Wxz*zdot*sin(chi)/(cos(gamma)*V^2));',f21,'f21',0);
    f22 = Wxz*cos(chi)*zdot/(V*cos(gamma));
    f22 = adigatorVarAnalyzer('f22 = Wxz*cos(chi)*zdot/(V*cos(gamma));',f22,'f22',0);
    f23 = 0.5*rho*S*Cl*V*sin(nu)*sin(gamma)/(m*(cos(gamma))^2) + Wxz*zdot*sin(chi)*sin(gamma)/(V*(cos(gamma))^2);
    f23 = adigatorVarAnalyzer('f23 = 0.5*rho*S*Cl*V*sin(nu)*sin(gamma)/(m*(cos(gamma))^2) + Wxz*zdot*sin(chi)*sin(gamma)/(V*(cos(gamma))^2);',f23,'f23',0);
    f24 = 0;
    f24 = adigatorVarAnalyzer('f24 = 0;',f24,'f24',0);
    f25 = 0;
    f25 = adigatorVarAnalyzer('f25 = 0;',f25,'f25',0);
    f26 = p_exp*(p_exp-1)*VR*((-z)^(p_exp-2))*zdot*sin(chi)/(V*cos(gamma));
    f26 = adigatorVarAnalyzer('f26 = p_exp*(p_exp-1)*VR*((-z)^(p_exp-2))*zdot*sin(chi)/(V*cos(gamma));',f26,'f26',0);
    f31 = 0.5*rho*S*Cl*cos(nu)/m + (g*cos(gamma)/V^2) - (Wxz*cos(chi)*sin(gamma)*zdot/V^2);
    f31 = adigatorVarAnalyzer('f31 = 0.5*rho*S*Cl*cos(nu)/m + (g*cos(gamma)/V^2) - (Wxz*cos(chi)*sin(gamma)*zdot/V^2);',f31,'f31',0);
    f32 = -Wxz*zdot*sin(chi)*sin(gamma)/V;
    f32 = adigatorVarAnalyzer('f32 = -Wxz*zdot*sin(chi)*sin(gamma)/V;',f32,'f32',0);
    f33 = g*sin(gamma)/V + (Wxz*zdot*cos(gamma)*cos(chi)/V);
    f33 = adigatorVarAnalyzer('f33 = g*sin(gamma)/V + (Wxz*zdot*cos(gamma)*cos(chi)/V);',f33,'f33',0);
    f34 = 0;
    f34 = adigatorVarAnalyzer('f34 = 0;',f34,'f34',0);
    f35 = 0;
    f35 = adigatorVarAnalyzer('f35 = 0;',f35,'f35',0);
    f36 = p_exp*(p_exp-1)*VR*((-z)^(p_exp-2))*cos(chi)*sin(gamma)*zdot/V;
    f36 = adigatorVarAnalyzer('f36 = p_exp*(p_exp-1)*VR*((-z)^(p_exp-2))*cos(chi)*sin(gamma)*zdot/V;',f36,'f36',0);
    f41 = cos(chi)*cos(gamma);
    f41 = adigatorVarAnalyzer('f41 = cos(chi)*cos(gamma);',f41,'f41',0);
    f42 = -V*sin(chi)*cos(gamma);
    f42 = adigatorVarAnalyzer('f42 = -V*sin(chi)*cos(gamma);',f42,'f42',0);
    f43 = -V*cos(chi)*sin(gamma);
    f43 = adigatorVarAnalyzer('f43 = -V*cos(chi)*sin(gamma);',f43,'f43',0);
    f44 = 0;
    f44 = adigatorVarAnalyzer('f44 = 0;',f44,'f44',0);
    f45 = 0;
    f45 = adigatorVarAnalyzer('f45 = 0;',f45,'f45',0);
    f46 = Wxz;
    f46 = adigatorVarAnalyzer('f46 = Wxz;',f46,'f46',0);
    f51 = sin(chi)*cos(gamma);
    f51 = adigatorVarAnalyzer('f51 = sin(chi)*cos(gamma);',f51,'f51',0);
    f52 = V*cos(chi)*cos(gamma);
    f52 = adigatorVarAnalyzer('f52 = V*cos(chi)*cos(gamma);',f52,'f52',0);
    f53 = -V*sin(chi)*sin(gamma);
    f53 = adigatorVarAnalyzer('f53 = -V*sin(chi)*sin(gamma);',f53,'f53',0);
    f54 = 0;
    f54 = adigatorVarAnalyzer('f54 = 0;',f54,'f54',0);
    f55 = 0;
    f55 = adigatorVarAnalyzer('f55 = 0;',f55,'f55',0);
    f56 = 0;
    f56 = adigatorVarAnalyzer('f56 = 0;',f56,'f56',0);
    f61 = -sin(gamma);
    f61 = adigatorVarAnalyzer('f61 = -sin(gamma);',f61,'f61',0);
    f62 = 0;
    f62 = adigatorVarAnalyzer('f62 = 0;',f62,'f62',0);
    f63 = -V*cos(gamma);
    f63 = adigatorVarAnalyzer('f63 = -V*cos(gamma);',f63,'f63',0);
    f64 = 0;
    f64 = adigatorVarAnalyzer('f64 = 0;',f64,'f64',0);
    f65 = 0;
    f65 = adigatorVarAnalyzer('f65 = 0;',f65,'f65',0);
    f66 = 0;
    f66 = adigatorVarAnalyzer('f66 = 0;',f66,'f66',0);
    A = [f11, f12, f13, f14, f15, f16;             f21, f22, f23, f24, f25, f26;             f31, f32, f33, f34, f35, f36;             f41, f42, f43, f44, f45, f46;             f51, f52, f53, f54, f55, f56;             f61, f62, f63, f64, f65, f66];
    A = adigatorVarAnalyzer('A = [f11, f12, f13, f14, f15, f16;             f21, f22, f23, f24, f25, f26;             f31, f32, f33, f34, f35, f36;             f41, f42, f43, f44, f45, f46;             f51, f52, f53, f54, f55, f56;             f61, f62, f63, f64, f65, f66];',A,'A',0);
adigatorIfIterEnd(1,1);
[adigatorIfEvalStr, adigatorIfEvalVar] = adigatorIfIterStart(1,2);%#ok<NASGU>
if ~isempty(adigatorIfEvalStr)
adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorIfEvalStr); adigatorSetCellEvalFlag(0);
end
    I = eye(6);
    I = adigatorVarAnalyzer('I = eye(6);',I,'I',0);
    A = zeros(6);
    A = adigatorVarAnalyzer('A = zeros(6);',A,'A',0);
    % ADiGator FOR Statement #2: START
    cadaforvar2 = 1:6;
    cadaforvar2 = adigatorVarAnalyzer('cadaforvar2 = 1:6;',cadaforvar2,'cadaforvar2',0);
    [adigatorForVariable2, adigatorForEvalStr, adigatorForEvalVar] = adigatorForInitialize(2,cadaforvar2,0);%#ok<NASGU>
    if ~isempty(adigatorForEvalStr)
        adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorForEvalStr); adigatorSetCellEvalFlag(0);
    end
    for adigatorForVariable2i = adigatorForVariable2
    cadaforcount2 = adigatorForIterStart(2,adigatorForVariable2i);
    i = cadaforvar2(:,cadaforcount2);
    i = adigatorVarAnalyzer('i = cadaforvar2(:,cadaforcount2);',i,'i',0);
        adigatorVarAnalyzer('% get states at current nominal point');
        % Call to User Function get_traj --- (FunID 5)
        cadainput5_1 = t;
        cadainput5_1 = adigatorVarAnalyzer('cadainput5_1 = t;',cadainput5_1,'cadainput5_1',0);
        cadainput5_2 = tf;
        cadainput5_2 = adigatorVarAnalyzer('cadainput5_2 = tf;',cadainput5_2,'cadainput5_2',0);
        cadainput5_3 = coeffs;
        cadainput5_3 = adigatorVarAnalyzer('cadainput5_3 = coeffs;',cadainput5_3,'cadainput5_3',0);
        cadainput5_4 = N;
        cadainput5_4 = adigatorVarAnalyzer('cadainput5_4 = N;',cadainput5_4,'cadainput5_4',0);
        adigatorInputs = {cadainput5_1;cadainput5_2;cadainput5_3;cadainput5_4};
        [adigatorFunInfo, adigatorOutputs] = adigatortempfunc5(adigatorFunInfo,adigatorInputs);
        cadaoutput5_1 = adigatorOutputs{1};
        sig = cadaoutput5_1;
        sig = adigatorVarAnalyzer('sig = cadaoutput5_1;',sig,'sig',0);
        % Call to User Function get_xu --- (FunID 3)
        cadainput3_1 = sig;
        cadainput3_1 = adigatorVarAnalyzer('cadainput3_1 = sig;',cadainput3_1,'cadainput3_1',0);
        cadainput3_2 = VR;
        cadainput3_2 = adigatorVarAnalyzer('cadainput3_2 = VR;',cadainput3_2,'cadainput3_2',0);
        adigatorInputs = {cadainput3_1;cadainput3_2};
        [adigatorFunInfo, adigatorOutputs] = adigatortempfunc3(adigatorFunInfo,adigatorInputs);
        cadaoutput3_1 = adigatorOutputs{1};
        X = cadaoutput3_1;
        X = adigatorVarAnalyzer('X = cadaoutput3_1;',X,'X',0);
        cadaoutput3_2 = adigatorOutputs{2};
        h = 1e-4;
        h = adigatorVarAnalyzer('h = 1e-4;',h,'h',0);
        x0 = [X(1); X(3); X(2); sig(1); sig(2); sig(3)];
        x0 = adigatorVarAnalyzer('x0 = [X(1); X(3); X(2); sig(1); sig(2); sig(3)];',x0,'x0',0);
        % Call to User Function non_flat_model --- (FunID 4)
        cadainput4_1 = t;
        cadainput4_1 = adigatorVarAnalyzer('cadainput4_1 = t;',cadainput4_1,'cadainput4_1',0);
        cadainput4_2 = x0 - h*I(:,i);
        cadainput4_2 = adigatorVarAnalyzer('cadainput4_2 = x0 - h*I(:,i);',cadainput4_2,'cadainput4_2',0);
        cadainput4_3 = solution;
        cadainput4_3 = adigatorVarAnalyzer('cadainput4_3 = solution;',cadainput4_3,'cadainput4_3',0);
        adigatorInputs = {cadainput4_1;cadainput4_2;cadainput4_3};
        [adigatorFunInfo, adigatorOutputs] = adigatortempfunc4(adigatorFunInfo,adigatorInputs);
        cadaoutput4_1 = adigatorOutputs{1};
        f1 = cadaoutput4_1;
        f1 = adigatorVarAnalyzer('f1 = cadaoutput4_1;',f1,'f1',0);
        % Call to User Function non_flat_model --- (FunID 4)
        cadainput4_1 = t;
        cadainput4_1 = adigatorVarAnalyzer('cadainput4_1 = t;',cadainput4_1,'cadainput4_1',0);
        cadainput4_2 = x0 + h*I(:,i);
        cadainput4_2 = adigatorVarAnalyzer('cadainput4_2 = x0 + h*I(:,i);',cadainput4_2,'cadainput4_2',0);
        cadainput4_3 = solution;
        cadainput4_3 = adigatorVarAnalyzer('cadainput4_3 = solution;',cadainput4_3,'cadainput4_3',0);
        adigatorInputs = {cadainput4_1;cadainput4_2;cadainput4_3};
        [adigatorFunInfo, adigatorOutputs] = adigatortempfunc4(adigatorFunInfo,adigatorInputs);
        cadaoutput4_1 = adigatorOutputs{1};
        f2 = cadaoutput4_1;
        f2 = adigatorVarAnalyzer('f2 = cadaoutput4_1;',f2,'f2',0);
        if ~exist('A','var'); A = cadastruct([],'A',[],0); end
        A(:,i) = (f2-f1)/(2*h);
        A = adigatorVarAnalyzer('A(:,i) = (f2-f1)/(2*h);',A,'A',1);
    [adigatorForEvalStr, adigatorForEvalVar]= adigatorForIterEnd(2,adigatorForVariable2i);
    if ~isempty(adigatorForEvalStr)
        adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorForEvalStr); adigatorSetCellEvalFlag(0);
    end
    end
    % ADiGator FOR Statement #2: END
[adigatorIfEvalStr, adigatorIfEvalVar] = adigatorIfIterEnd(1,2);%#ok<NASGU>
if ~isempty(adigatorIfEvalStr)
adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorIfEvalStr); adigatorSetCellEvalFlag(0);
end
% ADiGator IF Statement #1: END
adigatorOutputs = {A};
[adigatorFunInfo, adigatorOutputs] = adigatorFunctionEnd(2,adigatorFunInfo,adigatorOutputs);
