clear all
close all

load acmod
load('./results/loiter_lin_001')

ac.p_exp = 1;

N = (length(Z)-3)/8;

dsSoln = struct;
dsSoln.N = N;
dsSoln.p_exp = ac.p_exp;

VR = Z(end-1); dsSoln.VR = VR;
tf = Z(end); dsSoln.T = tf; dsSoln.t = linspace(tf/N,tf,N).';
chi0 = Z(end-2); dsSoln.chi0 = chi0;

V = Z((1-1)*N+1:1*N); dsSoln.V = V;
chi = Z((2-1)*N+1:2*N); dsSoln.chi = chi;
gam = Z((3-1)*N+1:3*N); dsSoln.gam = gam;
x = Z((4-1)*N+1:4*N); dsSoln.x = x; 
y = Z((5-1)*N+1:5*N); dsSoln.y = y;
z = Z((6-1)*N+1:6*N); dsSoln.z = z;
CL = Z((7-1)*N+1:7*N); dsSoln.CL = CL;
mu = Z((8-1)*N+1:8*N); dsSoln.mu = mu;


save('./results/loiterExp_001','dsSoln');



