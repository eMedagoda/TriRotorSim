% TriRotor simulation

clear all
% close all
% clc

global Params

TIME = 50;
DT = 0.01;
T_SIZE = TIME/DT;
RAD2DEG = (180.0/pi);
DEG2RAD = 1.0/RAD2DEG;

Params.g = 9.81;

Params.m = 3.00;   % mass
Params.h = 0.15; % height
Params.w = 0.15; % width
Params.d = 0.15; % depth

Params.L1 = 0.05; % moment arm length
Params.L2 = 0.21; % moment arm length
Params.L3 = 0.05; % moment arm length
Params.L4 = 0.21; % moment arm length
Params.L5 = 0.53; % moment arm length

Params.Ixx = (1/12) * Params.m * (Params.h * Params.h + Params.d * Params.d); % moment of inertia
Params.Iyy = (1/12) * Params.m * (Params.w * Params.w + Params.d * Params.d); % moment of inertia
Params.Izz = (1/12) * Params.m * (Params.w * Params.w + Params.h * Params.h); % moment of inertia
Params.Ixz = 0.00;

X = zeros(12,T_SIZE); % u,v,w,p,q,r,phi,theta,psi,Xdot,Ydot,Zdot,X,Y,Z
U = zeros(5,T_SIZE); % Fz, Fx, Mx, My, Mz
H = zeros(5,T_SIZE); % T_R, T_L, T_T, mu, d_mu

XTrim = [5.0,0.0,0.0*DEG2RAD]; % speed, alt, pitch

col = 'r';

% calculate trim states
[X0,U0,n] = Trim(XTrim(1),XTrim(2),XTrim(3))

% trim state rates
Xdot = TriSim_StateRates(X0,U0)

% trim motor states
H0 = TriSim_Motor(U0)

X(:,1) = X0;
U(1,:) = U0(1);
U(2,:) = U0(2);
U(3,:) = U0(3);
U(4,:) = U0(4);
U(5,:) = U0(5);
H(:,1) = H0;

[A,B] = TriSim_Linearise(X0,U0);
Alon = A([1,3,5,8,12],[1,3,5,8,12]);
Blon = B([1,3,5,8,12],[1,2,4]);
Alat = A([2,4,6,7,9],[2,4,6,7,9]);
Blat = B([2,4,6,7,9],[3,5]);

%##########################################################################

Clon = [1 0 0 0 0;0 0 0 0 -1];
Dlon = 0;
iu = [1, 2, 3];

tau_u      = 0.2;
tau_h      = 1.0;

PTv_all = [0 tau_u;0 tau_h]; % prediction horizon

nPH_all = size(PTv_all,2) - 1;  % no. of prediction steps

nOut_all = size(Clon,1);   % number of tracked output
nCont = length(iu);

% output weights (relative)
q_u      = 100.0;
q_h      = 100000.0;

% control weights (relative)
r_Fz = 10000.0;
r_Fx = 10000.0;
r_My = 10000.0;

q(1,1) = q_u;
q(2,1) = q_h;
Q = diag(repmat(q,[nPH_all,1]));        % final Q matrix (state output error)

r(1,1) = r_Fz;  
r(2,1) = r_Fx;
r(3,1) = r_My;
R = diag(repmat(r,[nPH_all,1]));        % final Q matrix (state output error)

MPCparams = ParametersMPCMulti(Alon,Blon,Clon,Q,R,iu,PTv_all);

MPCparams.K*MPCparams.F

Yc = [0.0;0.0];
Yref0(1,:) = repmat(Yc(1,1),1,MPCparams.nPH);
Yref0(2,:) = repmat(Yc(2,1),1,MPCparams.nPH);
Yref = reshape(Yref0,MPCparams.nPH*MPCparams.nOut,1);

%##########################################################################

for i = 2:T_SIZE
    
    X(:,i) = IntegrateRungeKutta(X(:,i-1),U(:,i-1),DT);  
    
    %######################################################################
    if (i >= 10/DT)        
       
        Yc = [1.0;0.0];
        
    end    
    
    Yref0(1,:) = repmat(Yc(1,1),1,MPCparams.nPH);
    Yref0(2,:) = repmat(Yc(2,1),1,MPCparams.nPH);
    Yref = reshape(Yref0,MPCparams.nPH*MPCparams.nOut,1);
    
    E = Yref - MPCparams.F * (X([1,3,5,8,12],i) - X0([1,3,5,8,12],1));
    
    UU = MPCparams.K * E;
    
    U(1,i) = U0(1) + UU(1); % Fz
    U(2,i) = U0(2) + UU(2); % Fx
    U(5,i) = U0(5) + UU(3); % My
    %######################################################################

%     dFz = 0.0;
%     
%     if (i >= 10/DT && i <= 11/DT)
%         dFz= -1.0;
%     end
%     
%     U(1,i) = U0(1) + dFz; % Fx
           
    H(:,i) = TriSim_Motor(U(:,i));    
    
end

t = DT:DT:TIME;

figure(1)
plot3(X(11,:),X(10,:),-X(12,:),col,'Linewidth',2)
xlim([-50 50])
ylim([-50 50])
zlim([-50 50])
grid
hold on

figure(2)
subplot(3,1,1)
plot(t,X(1,:),col,'Linewidth',2)
xlabel('Time (s)')
ylabel('u')
hold on
subplot(3,1,2)
plot(t,X(2,:),col,'Linewidth',2)
xlabel('Time (s)')
ylabel('v')
hold on
subplot(3,1,3)
plot(t,X(3,:),col,'Linewidth',2)
xlabel('Time (s)')
ylabel('w')
hold on

figure(3)
subplot(3,1,1)
plot(t,X(4,:)*RAD2DEG,col,'Linewidth',2)
xlabel('Time (s)')
ylabel('p')
hold on
subplot(3,1,2)
plot(t,X(5,:)*RAD2DEG,col,'Linewidth',2)
xlabel('Time (s)')
ylabel('q')
hold on
subplot(3,1,3)
plot(t,X(6,:)*RAD2DEG,col,'Linewidth',2)
xlabel('Time (s)')
ylabel('r')
hold on

figure(4)
subplot(3,1,1)
plot(t,X(7,:)*RAD2DEG,col,'Linewidth',2)
xlabel('Time (s)')
ylabel('phi')
hold on
subplot(3,1,2)
plot(t,X(8,:)*RAD2DEG,col,'Linewidth',2)
xlabel('Time (s)')
ylabel('theta')
hold on
subplot(3,1,3)
plot(t,X(9,:)*RAD2DEG,col,'Linewidth',2)
xlabel('Time (s)')
ylabel('psi')
hold on

figure(10)
subplot(3,1,1)
plot(t,X(10,:),col,'Linewidth',2)
xlabel('Time (s)')
ylabel('X')
hold on
subplot(3,1,2)
plot(t,X(11,:),col,'Linewidth',2)
xlabel('Time (s)')
ylabel('Y')
hold on
subplot(3,1,3)
plot(t,-X(12,:),col,'Linewidth',2)
xlabel('Time (s)')
ylabel('h')
hold on

figure(11)
subplot(3,1,1)
plot(t,H(1,:),col,'Linewidth',2)
xlabel('Time (s)')
ylabel('T_R')
hold on
subplot(3,1,2)
plot(t,H(2,:),col,'Linewidth',2)
xlabel('Time (s)')
ylabel('T_L')
hold on
subplot(3,1,3)
plot(t,H(3,:),col,'Linewidth',2)
xlabel('Time (s)')
ylabel('T_T')
hold on

figure(12)
subplot(2,1,1)
plot(t,H(4,:)*RAD2DEG,col,'Linewidth',2)
xlabel('Time (s)')
ylabel('rotor tilt')
hold on
subplot(2,1,2)
plot(t,H(5,:)*RAD2DEG,col,'Linewidth',2)
xlabel('Time (s)')
ylabel('rotor tilt offset')
hold on

figure(13)
subplot(5,1,1)
plot(t,U(1,:),col,'Linewidth',2)
hold on
subplot(5,1,2)
plot(t,U(2,:),col,'Linewidth',2)
hold on
subplot(5,1,3)
plot(t,U(3,:),col,'Linewidth',2)
hold on
subplot(5,1,4)
plot(t,U(4,:),col,'Linewidth',2)
hold on
subplot(5,1,5)
plot(t,U(5,:),col,'Linewidth',2)
hold on

figure(15)
plot3(X(11,:),X(10,:),-X(12,:),col,'Linewidth',2)
hold on
axis equal
