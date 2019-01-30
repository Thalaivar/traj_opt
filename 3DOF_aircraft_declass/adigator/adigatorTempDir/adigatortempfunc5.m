function [adigatorFunInfo, adigatorOutputs] = adigatortempfunc5(adigatorFunInfo,adigatorInputs)
[flag, adigatorFunInfo, adigatorInputs] = adigatorFunctionInitialize(5,adigatorFunInfo,adigatorInputs);
if flag; adigatorOutputs = adigatorInputs; return; end;
t = adigatorInputs{1};
tf = adigatorInputs{2};
coeffs = adigatorInputs{3};
N = adigatorInputs{4};
nargin = 4; nargout = 1;ax = coeffs(:,1);
ax = adigatorVarAnalyzer('ax = coeffs(:,1);',ax,'ax',0);
ay = coeffs(:,2);
ay = adigatorVarAnalyzer('ay = coeffs(:,2);',ay,'ay',0);
az = coeffs(:,3);
az = adigatorVarAnalyzer('az = coeffs(:,3);',az,'az',0);
x = ax(1);
x = adigatorVarAnalyzer('x = ax(1);',x,'x',0);
y = ay(1);
y = adigatorVarAnalyzer('y = ay(1);',y,'y',0);
z = az(1);
z = adigatorVarAnalyzer('z = az(1);',z,'z',0);
xdot = 0;
xdot = adigatorVarAnalyzer('xdot = 0;',xdot,'xdot',0);
ydot = 0;
ydot = adigatorVarAnalyzer('ydot = 0;',ydot,'ydot',0);
zdot = 0;
zdot = adigatorVarAnalyzer('zdot = 0;',zdot,'zdot',0);
xddot = 0;
xddot = adigatorVarAnalyzer('xddot = 0;',xddot,'xddot',0);
yddot = 0;
yddot = adigatorVarAnalyzer('yddot = 0;',yddot,'yddot',0);
zddot = 0;
zddot = adigatorVarAnalyzer('zddot = 0;',zddot,'zddot',0);
% ADiGator FOR Statement #2: START
cadaforvar2 = 2:N+1;
cadaforvar2 = adigatorVarAnalyzer('cadaforvar2 = 2:N+1;',cadaforvar2,'cadaforvar2',0);
[adigatorForVariable2, adigatorForEvalStr, adigatorForEvalVar] = adigatorForInitialize(2,cadaforvar2,0);%#ok<NASGU>
if ~isempty(adigatorForEvalStr)
    adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorForEvalStr); adigatorSetCellEvalFlag(0);
end
for adigatorForVariable2i = adigatorForVariable2
cadaforcount2 = adigatorForIterStart(2,adigatorForVariable2i);
i = cadaforvar2(:,cadaforcount2);
i = adigatorVarAnalyzer('i = cadaforvar2(:,cadaforcount2);',i,'i',0);
    x = x + ax(i)*cos(2*pi*(i-1)*t/tf) + ax(i+N)*sin(2*pi*(i-1)*t/tf);
    x = adigatorVarAnalyzer('x = x + ax(i)*cos(2*pi*(i-1)*t/tf) + ax(i+N)*sin(2*pi*(i-1)*t/tf);',x,'x',0);
    y = y + ay(i)*cos(2*pi*(i-1)*t/tf) + ay(i+N)*sin(2*pi*(i-1)*t/tf);
    y = adigatorVarAnalyzer('y = y + ay(i)*cos(2*pi*(i-1)*t/tf) + ay(i+N)*sin(2*pi*(i-1)*t/tf);',y,'y',0);
    z = z + az(i)*cos(2*pi*(i-1)*t/tf) + az(i+N)*sin(2*pi*(i-1)*t/tf);
    z = adigatorVarAnalyzer('z = z + az(i)*cos(2*pi*(i-1)*t/tf) + az(i+N)*sin(2*pi*(i-1)*t/tf);',z,'z',0);
    xdot = xdot - ax(i)*(2*pi*(i-1)/tf)*sin(2*pi*(i-1)*t/tf) + ax(i+N)*(2*pi*(i-1)/tf)*cos(2*pi*(i-1)*t/tf);
    xdot = adigatorVarAnalyzer('xdot = xdot - ax(i)*(2*pi*(i-1)/tf)*sin(2*pi*(i-1)*t/tf) + ax(i+N)*(2*pi*(i-1)/tf)*cos(2*pi*(i-1)*t/tf);',xdot,'xdot',0);
    ydot = ydot - ay(i)*(2*pi*(i-1)/tf)*sin(2*pi*(i-1)*t/tf) + ay(i+N)*(2*pi*(i-1)/tf)*cos(2*pi*(i-1)*t/tf);
    ydot = adigatorVarAnalyzer('ydot = ydot - ay(i)*(2*pi*(i-1)/tf)*sin(2*pi*(i-1)*t/tf) + ay(i+N)*(2*pi*(i-1)/tf)*cos(2*pi*(i-1)*t/tf);',ydot,'ydot',0);
    zdot = zdot - az(i)*(2*pi*(i-1)/tf)*sin(2*pi*(i-1)*t/tf) + az(i+N)*(2*pi*(i-1)/tf)*cos(2*pi*(i-1)*t/tf);
    zdot = adigatorVarAnalyzer('zdot = zdot - az(i)*(2*pi*(i-1)/tf)*sin(2*pi*(i-1)*t/tf) + az(i+N)*(2*pi*(i-1)/tf)*cos(2*pi*(i-1)*t/tf);',zdot,'zdot',0);
    xddot = xddot - ax(i)*((2*pi*(i-1)/tf)^2)*cos(2*pi*(i-1)*t/tf) - ax(i+N)*((2*pi*(i-1)/tf)^2)*sin(2*pi*(i-1)*t/tf);
    xddot = adigatorVarAnalyzer('xddot = xddot - ax(i)*((2*pi*(i-1)/tf)^2)*cos(2*pi*(i-1)*t/tf) - ax(i+N)*((2*pi*(i-1)/tf)^2)*sin(2*pi*(i-1)*t/tf);',xddot,'xddot',0);
    yddot = yddot - ay(i)*((2*pi*(i-1)/tf)^2)*cos(2*pi*(i-1)*t/tf) - ay(i+N)*((2*pi*(i-1)/tf)^2)*sin(2*pi*(i-1)*t/tf);
    yddot = adigatorVarAnalyzer('yddot = yddot - ay(i)*((2*pi*(i-1)/tf)^2)*cos(2*pi*(i-1)*t/tf) - ay(i+N)*((2*pi*(i-1)/tf)^2)*sin(2*pi*(i-1)*t/tf);',yddot,'yddot',0);
    zddot = zddot - az(i)*((2*pi*(i-1)/tf)^2)*cos(2*pi*(i-1)*t/tf) - az(i+N)*((2*pi*(i-1)/tf)^2)*sin(2*pi*(i-1)*t/tf);
    zddot = adigatorVarAnalyzer('zddot = zddot - az(i)*((2*pi*(i-1)/tf)^2)*cos(2*pi*(i-1)*t/tf) - az(i+N)*((2*pi*(i-1)/tf)^2)*sin(2*pi*(i-1)*t/tf);',zddot,'zddot',0);
[adigatorForEvalStr, adigatorForEvalVar]= adigatorForIterEnd(2,adigatorForVariable2i);
if ~isempty(adigatorForEvalStr)
    adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorForEvalStr); adigatorSetCellEvalFlag(0);
end
end
% ADiGator FOR Statement #2: END
sigma = [x,y,z,xdot,ydot,zdot,xddot,yddot,zddot];
sigma = adigatorVarAnalyzer('sigma = [x,y,z,xdot,ydot,zdot,xddot,yddot,zddot];',sigma,'sigma',0);
adigatorOutputs = {sigma};
[adigatorFunInfo, adigatorOutputs] = adigatorFunctionEnd(5,adigatorFunInfo,adigatorOutputs);
