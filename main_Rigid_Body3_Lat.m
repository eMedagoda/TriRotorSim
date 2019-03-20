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

theta_C = 0*DEG2RAD;
V_C = 0.01;
h_C = 0;

% calculate trim states
[X0,U0,n] = TrimRigidBody(V_C,h_C,theta_C);

% linearise system about current state
[A, B] = LineariseRigidBody(X0,U0);

%--------------------------------------------------------------------------

Acont = A([1,2,3,4,5,6,7,8,9,12],[1,2,3,4,5,6,7,8,9,12]);
Bcont = B([1,2,3,4,5,6,7,8,9,12],[1,2,3,4,5,6]);

Ccont = zeros(6,10);
Ccont(1,1) = 1; % u
Ccont(2,3) = 1; % w
Ccont(3,4) = 1; % p
Ccont(4,5) = 1; % q
Ccont(5,6) = 1; % r
Ccont(6,8) = 1; % theta

iu = [1, 3, 4, 5, 6];

tau_u = 0.2;
tau_w = 0.2;
tau_p = 1.0;
tau_q = 1.0;
tau_r = 1.0;
tau_theta  = 1.0;

PTv_all = [0 tau_u;
           0 tau_w;
           0 tau_p; 
           0 tau_q;
           0 tau_r;
           0 tau_theta]; % prediction horizon

nPH_all = size(PTv_all,2) - 1;  % no. of prediction steps

% output weights (relative)
q_u      = 100.0;
q_w      = 100.0;
q_p      = 10.0;
q_q      = 1.0;
q_r      = 10.0;
q_theta  = 1.0;

% control weights (relative)
r_Fx = 0.1;
r_Fz = 0.1;
r_Mx = 1;
r_My = 0.1;
r_Mz = 1;

q(1,1) = q_u;
q(2,1) = q_w;
q(3,1) = q_p;
q(4,1) = q_q;
q(5,1) = q_r;
q(6,1) = q_theta;
Q = diag(repmat(q,[nPH_all,1]));        % final Q matrix (state output error)

r(1,1) = r_Fx;  
r(2,1) = r_Fz;
r(3,1) = r_Mx;
r(4,1) = r_My;
r(5,1) = r_Mz;

sigma = 1;
R = sigma * diag(repmat(r,[nPH_all,1]));        % final Q matrix (state output error)

MPCparams = ParametersMPCMulti(Acont,Bcont,Ccont,Q,R,iu,PTv_all);

F0 = MPCparams.F;
% MPCparams.G
% MPCparams.K

%--------------------------------------------------------------------------

X(:,1) = X0;
U(1,:) = U0(1); % Fx
U(2,:) = U0(2); % Fy
U(3,:) = U0(3); % Fz
U(4,:) = U0(4); % Mx
U(5,:) = U0(5); % My
U(6,:) = U0(6); % Mz

C_bn = DirectionCosineMatrix(X(7,1),X(8,1),X(9,1));
Vs_C = 0.5*(h_C - X(12,1));
dV_C = C_bn * [0;0;Vs_C];

u_C = V_C*cos(X(8,1));
w_C = V_C*sin(X(8,1));

% inertial to body rate transformation
C_euler = [1  0          -sin(X(8,1));
           0  cos(X(7,1)) sin(X(7,1))*cos(X(8,1));
           0 -sin(X(7,1)) cos(X(7,1))*cos(X(8,1))];

rates = C_euler * [0.0*DEG2RAD;0.0*DEG2RAD;0.0*DEG2RAD];

Yc = [u_C + dV_C(1); w_C + dV_C(3); rates(1); rates(2); rates(3); theta_C];

% Yc = [V_C*cos(theta_C);V_C*sin(theta_C);0;0;0;theta_C;];
Yref0(1,:) = repmat(Yc(1,1),1,MPCparams.nPH);
Yref0(2,:) = repmat(Yc(2,1),1,MPCparams.nPH);
Yref0(3,:) = repmat(Yc(3,1),1,MPCparams.nPH);
Yref0(4,:) = repmat(Yc(4,1),1,MPCparams.nPH);
Yref0(5,:) = repmat(Yc(5,1),1,MPCparams.nPH);
Yref0(6,:) = repmat(Yc(6,1),1,MPCparams.nPH);
Yref = reshape(Yref0,MPCparams.nPH*MPCparams.nOut,1);

E = Yref - MPCparams.F * X([1,2,3,4,5,6,7,8,9,12],1);

UU = MPCparams.K * E;

U(1,1) = U0(1) + UU(1); % Fx
U(3,1) = U0(3) + UU(2); % Fz
U(4,1) = U0(4) + UU(3); % Mx
U(5,1) = U0(5) + UU(4); % My
U(6,1) = U0(6) + UU(5); % Mz

for i = 2:T_SIZE

    % update states
    [X(:,i), Xdot] = IntegrateRungeKuttaRigidBody(X(:,i-1),U(:,i-1),DT);
    X(9,i) = PiMinusPi(X(9,i));
    
    % current trim state
    V = norm(X([1:3],i));
    theta = X(8,i);
    h = X(12,i);
    
%     if (mod(i,(1/DT)) == 0)
        
        % re-linearise about current state
        [X0,U0,n] = TrimRigidBody(V,h,theta);
        [A, B] = LineariseRigidBody(X0,U0);
%         disp('lin')
    
%     end
    
    Acont = A([1,2,3,4,5,6,7,8,9,12],[1,2,3,4,5,6,7,8,9,12]);
    Bcont = B([1,2,3,4,5,6,7,8,9,12],[1,2,3,4,5,6]);
    
    % re-formulate controller
    MPCparams = ParametersMPCMulti(Acont,Bcont,Ccont,Q,R,iu,PTv_all);
    
    if (i >= 10/DT)
        
        theta_C = -15*DEG2RAD;
        V_C = 10.0;
        h_C = 5.0;
        
        C_bn = DirectionCosineMatrix(X(7,i),X(8,i),X(9,i));        
        Vs_C = 0.5*(h_C - X(12,i));        
        dV_C = C_bn * [0;0;Vs_C];
        
        u_C = V_C*cos(X(8,i));
        w_C = V_C*sin(X(8,i));
       
        % inertial to body rate transformation
        C_euler = [1  0          -sin(X(8,i));
                   0  cos(X(7,i)) sin(X(7,i))*cos(X(8,i));
                   0 -sin(X(7,i)) cos(X(7,i))*cos(X(8,i))];
                
        rates = C_euler * [0.0*DEG2RAD;0.0*DEG2RAD;0.0*DEG2RAD];
        
        Yc = [u_C + dV_C(1); w_C + dV_C(3); rates(1); rates(2); rates(3); theta_C];
        Yref0(1,:) = repmat(Yc(1,1),1,MPCparams.nPH);
        Yref0(2,:) = repmat(Yc(2,1),1,MPCparams.nPH);
        Yref0(3,:) = repmat(Yc(3,1),1,MPCparams.nPH);
        Yref0(4,:) = repmat(Yc(4,1),1,MPCparams.nPH);
        Yref0(5,:) = repmat(Yc(5,1),1,MPCparams.nPH);
        Yref0(6,:) = repmat(Yc(6,1),1,MPCparams.nPH);
        Yref = reshape(Yref0,MPCparams.nPH*MPCparams.nOut,1);
       
    end    
    
    E = Yref - MPCparams.F * X([1,2,3,4,5,6,7,8,9,12],i);
        
    UU = MPCparams.K * E;
    
    U(1,i) = U0(1) + UU(1); % Fx
    U(3,i) = U0(3) + UU(2); % Fz
    U(4,i) = U0(4) + UU(3); % Mx
    U(5,i) = U0(5) + UU(4); % My
    U(6,i) = U0(6) + UU(5); % Mz
    
%     U(1,i) = U0(1); % Fx
%     U(3,i) = U0(3); % Fz
%     U(4,i) = U0(4); % Mx
%     U(5,i) = U0(5); % My
%     U(6,i) = U0(6); % Mz
           
end

t = DT:DT:TIME;
col = 'b--';

figure(1)
plot3(X(11,:),X(10,:),X(12,:),col,'Linewidth',2)
grid
axis equal
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
