% TriRotor simulation

clear all
% close all
% clc

global Params

TIME = 30;
DT = 0.01;
T_SIZE = TIME/DT;
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
Params.Ixz = 0.0005;

X = zeros(12,T_SIZE); % u,v,w,p,q,r,phi,theta,psi,Xdot,Ydot,Zdot,X,Y,Z
U = zeros(6,T_SIZE);

theta_C = 0*pi/180;
V_C = 1.0;
h_C = 0;

% calculate trim states
[X0,U0,n] = TrimRigidBody(V_C,h_C,theta_C);

% linearise system about current state
[A, B] = LineariseRigidBody(X0,U0);

%--------------------------------------------------------------------------

Alat = A([2,4,6,7,9],[2,4,6,7,9]);
Blat = B([2,4,6,7,9],[2,4,6]);
Clat = [0 0 0 0 0;
        0 0 0 0 0;
        0 0 0 1 0;
        0 0 0 0 0];
    
Dlon = 0;
iu = [2, 3];

tau_v      = 0.5;
tau_phi    = 1.0;
tau_r      = 1.0;
tau_psi    = 1.0;

PTv_all = [0 tau_v;
           0 tau_phi;
           0 tau_r; 
           0 tau_psi]; % prediction horizon

nPH_all = size(PTv_all,2) - 1;  % no. of prediction steps

nOut_all = size(Clat,1);   % number of tracked output
nCont = length(iu);

% output weights (relative)
q_v      = 1.0;
q_phi    = 1.0;
q_r      = 10.0;
q_psi    = 1.0;

% control weights (relative)
r_Mx = 1.0;
r_Mz = 1.0;

q(1,1) = q_v;
q(2,1) = q_phi;
q(3,1) = q_r;
q(4,1) = q_psi;
Q = diag(repmat(q,[nPH_all,1]));        % final Q matrix (state output error)

r(1,1) = r_Mx;  
r(2,1) = r_Mz;
R = diag(repmat(r,[nPH_all,1]));        % final Q matrix (state output error)

MPCparams = ParametersMPCMulti(Alat,Blat,Clat,Q,R,iu,PTv_all);

% MPCparams.F
% MPCparams.G
% MPCparams.K

MPCparams.K*MPCparams.F

%--------------------------------------------------------------------------

Yc = [0.0;0.0;0.0;0.0];
Yref0(1,:) = repmat(Yc(1,1),1,MPCparams.nPH);
Yref0(2,:) = repmat(Yc(2,1),1,MPCparams.nPH);
Yref0(3,:) = repmat(Yc(3,1),1,MPCparams.nPH);
Yref0(4,:) = repmat(Yc(4,1),1,MPCparams.nPH);
Yref = reshape(Yref0,MPCparams.nPH*MPCparams.nOut,1);

X(:,1) = X0;
U(1,:) = U0(1); % Fx
U(2,:) = U0(2); % Fy
U(3,:) = U0(3); % Fz
U(4,:) = U0(4); % Mx
U(5,:) = U0(5); % My
U(6,:) = U0(6); % Mz

for i = 2:T_SIZE

    % update states
    [X(:,i), Xdot] = IntegrateRungeKuttaRigidBody(X(:,i-1),U(:,i-1),DT);
    X(9,i) = PiMinusPi(X(9,i));
    
    % current trim state
    V = norm(X([1:3],i));
    theta = X(8,i);
    h = X(12,i);
    
    % re-linearise about current state
    [X0,U0,n] = TrimRigidBody(V,h,theta);
    [A, B] = LineariseRigidBody(X0,U0);
    
    Alat = A([2,4,6,7,9],[2,4,6,7,9]);
    Blat = B([2,4,6,7,9],[2,4,6]);
    
    % re-formulate controller
    MPCparams = ParametersMPCMulti(Alat,Blat,Clat,Q,R,iu,PTv_all);
    
    r_C = 0.0;
    
    if (i >= 10/DT && i < 11/DT)        
       
        r_C = 10.0*DEG2RAD;        
        
    end    
    
    Yc = [0.0; 0.0;r_C; 0.0];
    Yref0(1,:) = repmat(Yc(1,1),1,MPCparams.nPH);
    Yref0(2,:) = repmat(Yc(2,1),1,MPCparams.nPH);
    Yref0(3,:) = repmat(Yc(3,1),1,MPCparams.nPH);
    Yref0(4,:) = repmat(Yc(4,1),1,MPCparams.nPH);
    Yref = reshape(Yref0,MPCparams.nPH*MPCparams.nOut,1);
    
    E = Yref - MPCparams.F * X([2,4,6,7,9],i);
    
    UU = MPCparams.K * E;
    
    U(4,i) = U0(4) + UU(1); % Mx
    U(6,i) = U0(6) + UU(2); % Mz
           
end

t = DT:DT:TIME;
col = 'r';

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

figure(5)
subplot(3,1,1)
plot(t,X(10,:),col,'Linewidth',2)
hold on
subplot(3,1,2)
plot(t,X(11,:),col,'Linewidth',2)
hold on
subplot(3,1,3)
plot(t,X(12,:),col,'Linewidth',2)
hold on

figure(13)
subplot(3,2,1)
plot(t,U(1,:),col,'Linewidth',2)
hold on
subplot(3,2,2)
plot(t,U(4,:),col,'Linewidth',2)
hold on
subplot(3,2,3)
plot(t,U(2,:),col,'Linewidth',2)
hold on
subplot(3,2,4)
plot(t,U(5,:),col,'Linewidth',2)
hold on
subplot(3,2,5)
plot(t,U(3,:),col,'Linewidth',2)
hold on
subplot(3,2,6)
plot(t,U(6,:),col,'Linewidth',2)
hold on
