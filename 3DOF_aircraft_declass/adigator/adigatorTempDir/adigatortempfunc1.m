function [adigatorFunInfo, adigatorOutputs] = adigatortempfunc1(adigatorFunInfo,adigatorInputs)
[flag, adigatorFunInfo, adigatorInputs] = adigatorFunctionInitialize(1,adigatorFunInfo,adigatorInputs);
if flag; adigatorOutputs = adigatorInputs; return; end;
X = adigatorInputs{1};
solution = adigatorInputs{2};
nargin = 2; nargout = 1;N = solution.N;
N = adigatorVarAnalyzer('N = solution.N;',N,'N',0);
n_coeffs = 2*N+1;
n_coeffs = adigatorVarAnalyzer('n_coeffs = 2*N+1;',n_coeffs,'n_coeffs',0);
adigatorVarAnalyzer('% to be used during trajectory optimisation');
% ADiGator IF Statement #1: START
cadaconditional1 = strcmp(solution.objfun_type, 'traj');
cadaconditional1 = adigatorVarAnalyzer('cadaconditional1 = strcmp(solution.objfun_type, ''traj'');',cadaconditional1,'cadaconditional1',0);
cadaconditional2 = strcmp(solution.objfun_type, 'floq_new');
cadaconditional2 = adigatorVarAnalyzer('cadaconditional2 = strcmp(solution.objfun_type, ''floq_new'');',cadaconditional2,'cadaconditional2',0);
adigatorIfInitialize(1,cadaconditional1,cadaconditional2);
adigatorIfIterStart(1,1);
    VR = X(3*n_coeffs+1,1);
    VR = adigatorVarAnalyzer('VR = X(3*n_coeffs+1,1);',VR,'VR',0);
    f = VR;
    f = adigatorVarAnalyzer('f = VR;',f,'f',0);
    adigatorVarAnalyzer('% to be used when trying new method of stability optimisation');
adigatorIfIterEnd(1,1);
[adigatorIfEvalStr, adigatorIfEvalVar] = adigatorIfIterStart(1,2);%#ok<NASGU>
if ~isempty(adigatorIfEvalStr)
adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorIfEvalStr); adigatorSetCellEvalFlag(0);
end
    coeffs_x = X(1:n_coeffs,1);
    coeffs_x = adigatorVarAnalyzer('coeffs_x = X(1:n_coeffs,1);',coeffs_x,'coeffs_x',0);
    coeffs_y = X(n_coeffs+1:2*n_coeffs,1);
    coeffs_y = adigatorVarAnalyzer('coeffs_y = X(n_coeffs+1:2*n_coeffs,1);',coeffs_y,'coeffs_y',0);
    coeffs_z = X(2*n_coeffs+1:3*n_coeffs,1);
    coeffs_z = adigatorVarAnalyzer('coeffs_z = X(2*n_coeffs+1:3*n_coeffs,1);',coeffs_z,'coeffs_z',0);
    coeffs = [coeffs_x,coeffs_y, coeffs_z];
    coeffs = adigatorVarAnalyzer('coeffs = [coeffs_x,coeffs_y, coeffs_z];',coeffs,'coeffs',0);
    tf = X(3*n_coeffs+2,1);
    tf = adigatorVarAnalyzer('tf = X(3*n_coeffs+2,1);',tf,'tf',0);
    VR = X(3*n_coeffs+1,1);
    VR = adigatorVarAnalyzer('VR = X(3*n_coeffs+1,1);',VR,'VR',0);
    if ~exist('solution','var'); solution = cadastruct([],'solution',[],0); end
    solution.tf = tf;
    solution = adigatorVarAnalyzer('solution.tf = tf;',solution,'solution',1);
    if ~exist('solution','var'); solution = cadastruct([],'solution',[],0); end
    solution.coeffs = coeffs;
    solution = adigatorVarAnalyzer('solution.coeffs = coeffs;',solution,'solution',1);
    if ~exist('solution','var'); solution = cadastruct([],'solution',[],0); end
    solution.VR = VR;
    solution = adigatorVarAnalyzer('solution.VR = VR;',solution,'solution',1);
    adigatorVarAnalyzer('% FTM by exponentials method');
    t = linspace(solution.tf, 0, 1000);
    t = adigatorVarAnalyzer('t = linspace(solution.tf, 0, 1000);',t,'t',0);
    FTM_expo = eye(3);
    FTM_expo = adigatorVarAnalyzer('FTM_expo = eye(3);',FTM_expo,'FTM_expo',0);
    del_t = t(1)-t(2);
    del_t = adigatorVarAnalyzer('del_t = t(1)-t(2);',del_t,'del_t',0);
    % ADiGator FOR Statement #1: START
    cadaforvar1 = 1:length(t);
    cadaforvar1 = adigatorVarAnalyzer('cadaforvar1 = 1:length(t);',cadaforvar1,'cadaforvar1',0);
    [adigatorForVariable1, adigatorForEvalStr, adigatorForEvalVar] = adigatorForInitialize(1,cadaforvar1,0);%#ok<NASGU>
    if ~isempty(adigatorForEvalStr)
        adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorForEvalStr); adigatorSetCellEvalFlag(0);
    end
    for adigatorForVariable1i = adigatorForVariable1
    cadaforcount1 = adigatorForIterStart(1,adigatorForVariable1i);
    i = cadaforvar1(:,cadaforcount1);
    i = adigatorVarAnalyzer('i = cadaforvar1(:,cadaforcount1);',i,'i',0);
        % Call to User Function get_jac --- (FunID 2)
        cadainput2_1 = solution;
        cadainput2_1 = adigatorVarAnalyzer('cadainput2_1 = solution;',cadainput2_1,'cadainput2_1',0);
        cadainput2_2 = t(i);
        cadainput2_2 = adigatorVarAnalyzer('cadainput2_2 = t(i);',cadainput2_2,'cadainput2_2',0);
        cadainput2_3 = 'FD';
        cadainput2_3 = adigatorVarAnalyzer('cadainput2_3 = ''FD'';',cadainput2_3,'cadainput2_3',0);
        adigatorInputs = {cadainput2_1;cadainput2_2;cadainput2_3};
        [adigatorFunInfo, adigatorOutputs] = adigatortempfunc2(adigatorFunInfo,adigatorInputs);
        cadaoutput2_1 = adigatorOutputs{1};
        Jk = cadaoutput2_1;
        Jk = adigatorVarAnalyzer('Jk = cadaoutput2_1;',Jk,'Jk',0);
        Jk = Jk(1:3,1:3);
        Jk = adigatorVarAnalyzer('Jk = Jk(1:3,1:3);',Jk,'Jk',0);
        FTM_expo = FTM_expo*expm(Jk*del_t);
        FTM_expo = adigatorVarAnalyzer('FTM_expo = FTM_expo*expm(Jk*del_t);',FTM_expo,'FTM_expo',0);
    [adigatorForEvalStr, adigatorForEvalVar]= adigatorForIterEnd(1,adigatorForVariable1i);
    if ~isempty(adigatorForEvalStr)
        adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorForEvalStr); adigatorSetCellEvalFlag(0);
    end
    end
    % ADiGator FOR Statement #1: END
    D = eig(FTM_expo);
    D = adigatorVarAnalyzer('D = eig(FTM_expo);',D,'D',0);
    f = 0;
    f = adigatorVarAnalyzer('f = 0;',f,'f',0);
    param1 = 1;
    param1 = adigatorVarAnalyzer('param1 = 1;',param1,'param1',0);
    param2 = 10;
    param2 = adigatorVarAnalyzer('param2 = 10;',param2,'param2',0);
    % ADiGator FOR Statement #2: START
    cadaforvar2 = 1:3;
    cadaforvar2 = adigatorVarAnalyzer('cadaforvar2 = 1:3;',cadaforvar2,'cadaforvar2',0);
    [adigatorForVariable2, adigatorForEvalStr, adigatorForEvalVar] = adigatorForInitialize(2,cadaforvar2,0);%#ok<NASGU>
    if ~isempty(adigatorForEvalStr)
        adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorForEvalStr); adigatorSetCellEvalFlag(0);
    end
    for adigatorForVariable2i = adigatorForVariable2
    cadaforcount2 = adigatorForIterStart(2,adigatorForVariable2i);
    i = cadaforvar2(:,cadaforcount2);
    i = adigatorVarAnalyzer('i = cadaforvar2(:,cadaforcount2);',i,'i',0);
        f = f + atan(param1*(abs(D(i)) - 1));
        f = adigatorVarAnalyzer('f = f + atan(param1*(abs(D(i)) - 1));',f,'f',0);
    [adigatorForEvalStr, adigatorForEvalVar]= adigatorForIterEnd(2,adigatorForVariable2i);
    if ~isempty(adigatorForEvalStr)
        adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorForEvalStr); adigatorSetCellEvalFlag(0);
    end
    end
    % ADiGator FOR Statement #2: END
    f = f + param2*VR;
    f = adigatorVarAnalyzer('f = f + param2*VR;',f,'f',0);
[adigatorIfEvalStr, adigatorIfEvalVar] = adigatorIfIterEnd(1,2);%#ok<NASGU>
if ~isempty(adigatorIfEvalStr)
adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorIfEvalStr); adigatorSetCellEvalFlag(0);
end
% ADiGator IF Statement #1: END
adigatorVarAnalyzer('%     % to be used when estimating floquet expo via collocation (probably');
adigatorVarAnalyzer('%     % obsolete)');
adigatorVarAnalyzer('%     elseif nargin > 3 && strcmp(type, ''floq_old'')');
adigatorVarAnalyzer('%         n = 3; % dimension of 3DOF system');
adigatorVarAnalyzer('%         eigval = complex(X(2*n*M+1,1), X(2*n*M+2,1));');
adigatorVarAnalyzer('%         f = -real(eigval);');
adigatorVarAnalyzer('%     end');
adigatorOutputs = {f};
[adigatorFunInfo, adigatorOutputs] = adigatorFunctionEnd(1,adigatorFunInfo,adigatorOutputs);
