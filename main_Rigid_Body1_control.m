% TriRotor simulation

% close all
% clc

global Params

TIME = 100;
DT = 0.1;
T_SIZE = TIME/DT;
RAD2DEG = (180/pi);
DEG2RAD = 1/RAD2DEG;

Params.g = 9.81;

Params.m = 3.0;   % mass
Params.h = 0.15; % height
Params.w = 0.15; % width
Params.d = 0.3; % depth

Params.Ax = Params.h*Params.w;
Params.Ay = Params.h*Params.d;
Params.Az = Params.w*Params.d;

Params.L1 = 0.05; % moment arm length
Params.L2 = 0.21; % moment arm length
Params.L3 = 0.05; % moment arm length
Params.L4 = 0.21; % moment arm length
Params.L5 = 0.53; % moment arm length

Params.Ixx = (1/12) * Params.m * (Params.h * Params.h + Params.d * Params.d); % moment of inertia
Params.Iyy = (1/12) * Params.m * (Params.w * Params.w + Params.d * Params.d); % moment of inertia
Params.Izz = (1/12) * Params.m * (Params.w * Params.w + Params.h * Params.h); % moment of inertia
Params.Ixz = 0.001;

X = zeros(12,T_SIZE); % u,v,w,p,q,r,phi,theta,psi,Xdot,Ydot,Zdot,X,Y,Z
U = zeros(6,T_SIZE);

V_C = 0.01;
h_C = 10.0;
theta_C = 0.0*DEG2RAD;
psi_C = 0.0*DEG2RAD;

% calculate trim states
[X0,U0,n] = TrimRigidBody(V_C,-h_C,theta_C);

% linearise system about current state
[A, B] = LineariseRigidBody(X0,U0);

X(:,1) = X0;
U(1,:) = U0(1); % Fx
U(2,:) = U0(2); % Fy
U(3,:) = U0(3); % Fz
U(4,:) = U0(4); % Mx
U(5,:) = U0(5); % My
U(6,:) = U0(6); % Mz

U_MOTOR = [U0(3) U0(1) U0(4) U0(5) U0(6)]';

H = zeros(5,T_SIZE);
H(:,1) = TriSim_Motor(U_MOTOR);

u_sum = 0.0;
w_sum = 0.0;

for i = 2:T_SIZE

    % update states
    [X(:,i), Xdot] = IntegrateRungeKuttaRigidBody(X(:,i-1),U(:,i-1),DT);
 
    X(9,i) = PiMinusPi(X(9,i)); 
    
    % ---------------------------- GUIDANCE ------------------------------------
    
    V = norm(X([1:3],i));
    
    dFx = 0.0;
    dFz = 0.0;
    dMx = 0.0;
    dMy = 0.0;
    dMz = 0.0;    
   
    if (i > 10/DT)
        
        V_C = 0.01;
        h_C = 10.0;
        theta_C = -45.0*DEG2RAD;
        psi_C = 0.0*pi/180;
       
    end
    
     if (i > 50/DT)
         h_C = 10.0;
     end
    
    % speed demand generation (navigation frame to body frame)
    C_bn = DirectionCosineMatrix(X(7,i),X(8,i),X(9,i));
    Vz_C = 1.0*(-h_C - X(12,i));
    dV_C = C_bn * [V_C*cos(X(9,i));V_C*sin(X(9,i));Vz_C];
    u_C = dV_C(1);
    v_C = dV_C(2);
    w_C = dV_C(3);
    
    %--------------------------- CONTROLLER ------------------------------------
    
    % forward speed control
    int_u = 0.07;
    e_u = u_C - X(1,i);
    dFx = 3.0 * e_u + int_u*u_sum; % apply intergral action on speed controllers
    
    % vertical speed control
    int_w = 0.07;
    e_w = w_C - X(3,i);
    dFz = 3.0 * e_w + int_w*w_sum;
    
    % pitch control
    e_theta = theta_C - X(8,i);
    dMy = 0.003 * e_theta + 0.01 * (0.0 - X(5,i));
    
    % roll regulator
    phi_C = 0.0*pi/180;
    e_phi = phi_C - X(7,i);
    dMx = 0.01 * e_phi + 0.02 * (0.0 - X(4,i)) + 0.0003 * (v_C - X(2,i)); % use bank to manage lateral velocity (side-slip)
    
    % heading controller
    e_psi = PiMinusPi(psi_C - X(9,i));    
    dMz = 0.01 * e_psi + 0.05 * (0.0 - X(6,i));
        
    % speed integrator engage
    if (V >= 0.95 * V_C)
        
        u_sum = u_sum + e_u;
        w_sum = w_sum + e_w;
        
    end
    
    %----------------- APPLY CONTROL OUTPUTS TO TRIM ---------------------------
   
    U(1,i) = U0(1) + dFx; % Fx
    U(3,i) = U0(3) + dFz; % Fz
    U(4,i) = U0(4) + dMx; % Mx
    U(5,i) = U0(5) + dMy; % My
    U(6,i) = U0(6) + dMz; % Mz
    
    U_MOTOR = [U(3,i) U(1,i) U(4,i) U(5,i) U(6,i)]';
    
    H(:,i) = TriSim_Motor(U_MOTOR);
    
end

t = DT:DT:TIME;
col = 'b';

figure(1)
plot3(X(11,:),X(10,:),-X(12,:),col,'Linewidth',2)
grid
axis equal
hold on

figure(2)
subplot(4,1,1)
plot(t,X(1,:),col,'Linewidth',2)
xlabel('Time (s)')
ylabel('u')
hold on
subplot(4,1,2)
plot(t,X(2,:),col,'Linewidth',2)
xlabel('Time (s)')
ylabel('v')
hold on
subplot(4,1,3)
plot(t,X(3,:),col,'Linewidth',2)
xlabel('Time (s)')
ylabel('w')
hold on
subplot(4,1,4)
plot(t,sqrt(X(1,:).*X(1,:) + X(2,:).*X(2,:) + X(3,:).*X(3,:)),col,'Linewidth',2)
xlabel('Time (s)')
ylabel('V')
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
xlabel('Time (s)')
ylabel('X(m)')
hold on
subplot(3,1,2)
plot(t,X(11,:),col,'Linewidth',2)
xlabel('Time (s)')
ylabel('Y(m)')
hold on
subplot(3,1,3)
plot(t,-X(12,:),col,'Linewidth',2)
xlabel('Time (s)')
ylabel('Alt(m)')
hold on

figure(6)
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

figure(7)
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

figure(8)
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

figure(9)
plot(X(11,:),X(10,:),col,'Linewidth',2)
grid
axis equal
hold on
