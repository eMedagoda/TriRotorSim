% Quadcopter simulation

clear all
close all
clc

global Params

RAD2DEG = (180/pi);
DEG2RAD = 1/RAD2DEG;

Params.g = 9.81;

Params.m = 1;   % mass
Params.h = 0.1; % height
Params.w = 0.1; % width
Params.d = 0.1; % depth

Params.L1 = 0.05; % moment arm length
Params.L2 = 0.21; % moment arm length
Params.L3 = 0.05; % moment arm length
Params.L4 = 0.21; % moment arm length
Params.L5 = 0.53; % moment arm length

Params.Ixx = (1/12) * Params.m * (Params.h * Params.h + Params.d * Params.d); % moment of inertia
Params.Iyy = (1/12) * Params.m * (Params.w * Params.w + Params.d * Params.d); % moment of inertia
Params.Izz = (1/12) * Params.m * (Params.w * Params.w + Params.h * Params.h); % moment of inertia
Params.Ixz = 0.001;

XTrim = [0,0,0*DEG2RAD];

col = 'b';

[X0,U0,n] = Trim(XTrim(1),XTrim(2),XTrim(3));
Xdot = TriSim_StateRates(X0,U0);

[A,B] = TriSim_Linearise(X0,U0);

ALon = A([1,3,5,8,12],[1,3,5,8,12]);
BLon = B([1,3,5,8,12],[1,2,4]);

ALat = A([2,4,6,7,9],[2,4,6,7,9]);
BLat = B([2,4,6,7,9],[3,5]);

[z,p,k] = ss2zp(ALon,BLon,[0 0 0 0 -1],[0 0 0],1);
G_U1_Z = minreal(zpk(z,p,k));

[z,p,k] = ss2zp(ALon,BLon,[0 1 0 0 0],[0 0 0],1);
G_U1_w = minreal(zpk(z,p,k));

[z,p,k] = ss2zp(ALon,BLon,[1 0 0 0 0],[0 0 0],2);
G_U2_u = minreal(zpk(z,p,k));

[z,p,k] = ss2zp(ALat,BLat,[0 1 0 0 0],[0 0],1);
G_U3_p = minreal(zpk(z,p,k));

[z,p,k] = ss2zp(ALat,BLat,[0 0 0 1 0],[0 0],1);
G_U3_phi = minreal(zpk(z,p,k));

[z,p,k] = ss2zp(ALat,BLat,[0 0 1 0 0],[0 0],2);
G_U5_r = minreal(zpk(z,p,k));

[z,p,k] = ss2zp(ALat,BLat,[0 0 0 0 1],[0 0],2);
G_U5_psi = minreal(zpk(z,p,k));

s = tf('s');
K_p = 0.01;
K_phi = 10;

G_p_phi = minreal(G_U3_phi/G_U3_p);
G_pC_p = minreal((G_U3_p)/(1 + (K_p*G_U3_p)));

G_pC_phi = minreal(G_p_phi*G_pC_p);

G_phiC_phi = minreal((K_phi*G_pC_phi)/(1 + (K_phi*G_pC_phi)))

step(G_phiC_phi)
