% TriRotor simulation

clear all
% close all
% clc

global Params

TIME = 80;
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
Params.Ixz = 0.001;

X = zeros(12,T_SIZE); % u,v,w,p,q,r,phi,theta,psi,Xdot,Ydot,Zdot,X,Y,Z
U = zeros(5,T_SIZE);

X_v = zeros(3,T_SIZE);
X_v(1,1) = 1;
X_v(2,1) = 0;
X_v(3,1) = 45*DEG2RAD;

XTrim = [1,0,0*DEG2RAD];

col = 'r';
col2 = 'k';

[X0,U0,n] = Trim(XTrim(1),XTrim(2),XTrim(3));
Xdot = TriSim_StateRates(X0,U0);

[X1,U1,n] = Trim(1,0,-20*DEG2RAD);
H1 = TriSim_Motor(U1);

X(:,1) = X0;
U(1,:) = U0(1);
U(2,:) = U0(2);
U(3,:) = U0(3);
U(4,:) = U0(4);
U(5,:) = U0(5);

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
K_p = 0.02;
K_phi = 0.8;
K_q = 0.2;
K_theta = 0.8;
K_r = 0.02;
K_psi = 0.8;

% K_u = zpk(50 * (s + 2)/((s + 9)*(s + 0.3)));
K_u = zpk(20 /(s + 5));
[A_u, B_u, C_u, D_u] = ssdata(K_u);
X_u = [0];

% K_w = zpk(50 * (s + 2)/((s + 9)*(s + 0.3)));
K_w = zpk(100 /(s + 5));
[A_w, B_w, C_w, D_w] = ssdata(K_w);
X_w = [0];

% K_h = zpk(-20 * (s + 1)/(s + 6));
K_h = zpk(-20 * (s + 5)/(s + 0.1));
[A_h, B_h, C_h, D_h] = ssdata(K_h);
X_h = 0;

%G_p_phi = minreal(G_U3_phi/G_U3_p);
%G_pC_p = minreal((K_p*G_U3_p)/(1 + (K_p*G_U3_p)));
%
%G_pC_phi = minreal(G_p_phi*G_pC_p);
%
%G_phiC_phi = minreal((K_phi*G_pC_phi)/(1 + (K_phi*G_pC_phi)));

H(:,1) = TriSim_Motor(U(:,1));

phi_C = 0*DEG2RAD;
theta_C = 0*DEG2RAD;
psi_C = 0*DEG2RAD;
h_C = 0;
V_C = 0;
u_C = V_C * cos(theta_C);
w_C = V_C * sin(theta_C);

Vg(1) = sqrt(Xdot(10)*Xdot(10) + Xdot(11)*Xdot(11));
Psi_dot(1) = Xdot(9);

e_u_sum = 0;
e_w_sum = 0;
e_h_sum = 0;

P_V = [X(10,1);X(11,1)];
P_T = [X_v(1,1);X_v(2,1)];
[V_C(1), psi_C(1), R(1)] = TriSim_Guidance(P_V,P_T);

for i = 2:T_SIZE

    
    [X(:,i), Xdot] = IntegrateRungeKutta(X(:,i-1),U(:,i-1),DT);
    X(9,i) = PiMinusPi(X(9,i));
    
    X_v(:,i) = VehSim_IntegrateRungeKutta(X_v(:,i-1),[1;30*DEG2RAD],DT);
    X_v(3,i) = PiMinusPi(X_v(3,i));
        
    Vg(i) = sqrt(Xdot(10)*Xdot(10) + Xdot(11)*Xdot(11));
    Psi_dot(i) = Xdot(9);
    
    %--------------------------------------------------
%     if i > 2/DT && i <= 15/DT
%         
%         V_C = 0;
%         theta_C = -60*DEG2RAD;
%         u_C = V_C * cos(theta_C);
%         w_C = V_C * sin(theta_C);
%         h_C = 0;
%         
%     elseif i > 15/DT
%         
%         V_C = 1;
%         theta_C = -60*DEG2RAD;
%         u_C = V_C * cos(theta_C);
%         w_C = V_C * sin(theta_C);        
%         
%     end
%     
%     if i > 4/DT && i <= 20/DT
%         
%         psi_C = 0*DEG2RAD;
%         
%     elseif i > 20/DT
%         
%         psi_C = -60*DEG2RAD;
%          
%     end

    if i > 15/DT
        
        theta_C = 0*DEG2RAD;
               
    end
    
    P_V = [X(10,i);X(11,i)];
%     P_T = [X_v(1,i);X_v(2,i)];    
    P_T = [10;10];

%     if i > 4/DT        
    
        [V_C(i), psi_C(i), R(i)] = TriSim_Guidance(P_V,P_T);
%         [V_C(i), psi_C(i)] = TriSim_Guidance2(P_V,P_T,0*DEG2RAD);
        u_C = V_C(i) * cos(theta_C);
        w_C = V_C(i) * sin(theta_C);
        
%     end
    
    %--------------------------------------------------
    
    X_h_dot = A_h * X_h + B_h * (h_C + X(12,i));
    X_h = X_h + X_h_dot*DT;
    e_h_sum = e_u_sum + (h_C + X(12,i));
    UU1h = C_h * X_h + D_h * (h_C + X(12,i)) + 0.0*e_h_sum;    
    U(1,i) = U(1,i) + UU1h; 
    
    X_w_dot = A_w * X_w + B_w * (w_C - X(3,i));
    X_w = X_w + X_w_dot*DT;
    e_w_sum = e_w_sum + (w_C - X(3,i));
    UU1w = C_w * X_w + D_w * (w_C - X(3,i)) + 0.05*e_w_sum;    
    U(1,i) = U(1,i) + UU1w;
    
    X_u_dot = A_u * X_u + B_u * (u_C - X(1,i));
    X_u = X_u + X_u_dot*DT;
    e_u_sum = e_u_sum + (u_C - X(1,i));
    UU2 = C_u * X_u + D_u * (u_C - X(1,i)) + 0.05*e_u_sum;    
    U(2,i) = U(2,i) + UU2;
    
    pC = K_phi*(phi_C - X(7,i));
    UU3 = K_p*(pC - X(4,i));
    U(3,i) = U(3,i) + UU3;
    
    qC = K_phi*(theta_C - X(8,i));
    UU4 = K_p*(qC - X(5,i));
    U(4,i) = U(4,i) + UU4;
    
    rC = K_phi*PiMinusPi(psi_C(i) - X(9,i));
    UU5 = K_p*(rC - X(6,i));
    U(5,i) = U(5,i) + UU5;
         
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
plot(t,Psi_dot*RAD2DEG,'k')

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
plot(t,psi_C*RAD2DEG,col2,'Linewidth',2)

figure(10)
subplot(3,1,1)
plot(t,X(10,:),col,'Linewidth',2)
hold on
subplot(3,1,2)
plot(t,X(11,:),col,'Linewidth',2)
hold on
subplot(3,1,3)
plot(t,-X(12,:),col,'Linewidth',2)
hold on

figure(11)
subplot(2,2,1)
plot(t,H(1,:),col,'Linewidth',2)
hold on
subplot(2,2,2)
plot(t,H(2,:),col,'Linewidth',2)
hold on
subplot(2,2,3)
plot(t,H(2,:)-H(1,:),col,'Linewidth',2)
hold on
subplot(2,2,4)
plot(t,H(3,:),col,'Linewidth',2)
hold on

figure(12)
subplot(2,1,1)
plot(t,H(4,:)*RAD2DEG,col,'Linewidth',2)
hold on
subplot(2,1,2)
plot(t,H(5,:)*RAD2DEG,col,'Linewidth',2)
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

V = sqrt(X(1,:).^2 + X(2,:).^2 + X(3,:).^2);
figure(14)
plot(t,V,col,'Linewidth',2)
hold on
plot(t,Vg,'k','Linewidth',2)
plot(t,V_C,'m','Linewidth',2)

figure(15)
plot(X(11,:),X(10,:),col,'Linewidth',2)
hold on
plot(X_v(2,:),X_v(1,:),col2,'Linewidth',2)
axis equal

figure(16)
plot(R,col,'Linewidth',2)
hold on
