% 30/05/19
close all
clear all

load acmod
load('./results/loiterExp_003');

ac.p_exp = 0.25;

N = (length(Z)-3)/18;

dsSoln = struct;
dsSoln.N = N;
dsSoln.p_exp = ac.p_exp;

VR = Z(end-1); dsSoln.VR = VR;
tf = Z(end); dsSoln.T = tf; dsSoln.t = linspace(tf/N,tf,N).';
Psi0 = Z(end-2); dsSoln.Psi0 = Psi0;

% X = [u,v,w,p,q,r,Phi,Thet,Psi,x,y,z,df,da,de,dr,CTx,CTy]

u =    Z((1-1)*N+1:1*N);   dsSoln.u = u;
v =    Z((2-1)*N+1:2*N);   dsSoln.v = v;
w =    Z((3-1)*N+1:3*N);   dsSoln.w = w;
p =    Z((4-1)*N+1:4*N);   dsSoln.p = p;
q =    Z((5-1)*N+1:5*N);   dsSoln.q = q;
r =    Z((6-1)*N+1:6*N);   dsSoln.r = r;
Phi =  Z((7-1)*N+1:7*N);   dsSoln.Phi = Phi;
Thet = Z((8-1)*N+1:8*N);   dsSoln.Thet = Thet;
Psi1 = Z((9-1)*N+1:9*N);   dsSoln.Psi1 = Psi1;
x =    Z((10-1)*N+1:10*N); dsSoln.x = x;
y =    Z((11-1)*N+1:11*N); dsSoln.y = y;
z =    Z((12-1)*N+1:12*N); dsSoln.z = z;

df =    Z((13-1)*N+1:13*N); dsSoln.df = df;
da =    Z((14-1)*N+1:14*N); dsSoln.da = da;
de =    Z((15-1)*N+1:15*N); dsSoln.de = de;
dr =    Z((16-1)*N+1:16*N); dsSoln.dr = dr;
CTx =    Z((17-1)*N+1:17*N); dsSoln.CTx = CTx;
CTy =    Z((18-1)*N+1:18*N); dsSoln.CTy = CTy;

save('./results/loiterExp_struct_003','dsSoln');



