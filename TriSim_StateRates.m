function Xdot = TriSim_StateRates(X,U)

global Params

g = Params.g;

Ixx = Params.Ixx;
Iyy = Params.Iyy;
Izz = Params.Izz;
Ixz = Params.Ixz;

m = Params.m;

% calculate body axis rotation rate coefficients
c0=Ixx*Izz-Ixz^2;
c1=Izz/c0;
c2=Ixz/c0;
c3=c2*(Ixx-Iyy+Izz);
c4=c1*(Iyy-Izz)-c2*Ixz;
c5=1/Iyy;
c6=c5*Ixz;
c7=c5*(Izz-Ixx);
c8=Ixx/c0;
c9=c8*(Ixx-Iyy)+c2*Ixz;

Xdot = zeros(12,1);

u = X(1);
v = X(2);
w = X(3);

p = X(4);
q = X(5);
r = X(6);

phi = X(7);
theta = X(8);
psi = X(9);

c = 0.1;
b = 0.1;
S = Params.d * Params.w;
Vt = sqrt(u^2 + v^2 + w^2);
rho = Atmosphere(-X(12,1));
qbar   = 0.5*rho*Vt^2;
bo2V   = 0.5*b/Vt;
co2V   = 0.5*c/Vt;

if u > 0.1
    alpha = atan(w/u);
    beta = atan(v/u);
    p_a  = p*bo2V;
    q_a  = q*co2V;
    r_a  = r*bo2V;
else
    alpha = 0;
    beta = 0;
    p_a = 0;
    q_a = 0;
    r_a = 0;
end

% alpha_o = -0.0*pi/180; 
% Cdo    =  1.500;
% k      =  0.050;
% CLa  =  5.827;
% CLq  =  7.960;
% CLo  = -CLa*alpha_o;
% Cmo  =  0.06;
% Cma  = -0.802;
% Cmq  = -17.72;
% Cyb  = -0.507;
% Cyp  = -0.128;
% Cyr  =  0.336;
% Cnb  =  0.107;
% Cnp  = -0.0226;
% Cnr  = -0.160;
% Clb  = -0.0852;
% Clp  = -0.328;
% Clr  =  0.0776;
% CL = CLo + CLa*alpha;
% Cd = Cdo + k*CL^2;

alpha_o = -0.0*pi/180; 
Cdo    =  1.05;
k      =  0.050;
CLa  =  5.827;
CLq  =  7.960;
CLo  = -CLa*alpha_o;
Cmo  =  0.06;
Cma  = -0.802;
Cmq  = -17.72;
Cyb  = -0.507;
Cyp  = -0.128;
Cyr  =  0.336;
Cnb  =  0.107;
Cnp  = -0.0226;
Cnr  = -0.160;
Clb  = -0.0852;
Clp  = -0.328;
Clr  =  0.0776;
CL = CLo + CLa*alpha;
Cd = Cdo + k*CL^2;

% force coefficients
Cx = -Cd*cos(alpha) + CL*sin(alpha);
Cy =  Cyb*beta + Cyp*p_a + Cyr*r_a;
Cz = -CL*cos(alpha) - Cd*sin(alpha) - CLq*q_a;

% moment coefficients
Cl = Clb*beta + Clp*p_a + Clr*r_a;
Cm = Cmo + Cma*alpha + Cmq*q_a;
Cn = Cnb*beta + Cnp*p_a + Cnr*r_a;

F_x = U(2) + qbar*S*Cx;
F_y = qbar*S*Cy;
F_z = U(1) + qbar*S*Cz;

M_x = U(3) + qbar*S*b*Cl;
M_y = U(4) + qbar*S*c*Cm;
M_z = U(5) + qbar*S*b*Cn;

% F_x = U(2) - sign(u)*0.5*u*u;
% F_y = - sign(v)*50*v*v;
% F_z = U(1) - sign(w)*1*w*w;
% 
% M_x = U(3);
% M_y = U(4);
% M_z = U(5);

udot = -g*sin(theta)            + r*v - q*w + F_x/m ;
vdot =  g*sin(phi)*cos(theta)   - r*u + p*w + F_y/m ;
wdot =  g*cos(phi)*cos(theta)   + q*u - p*v + F_z/m ;

% calculate body axis rotation rates
pdot = c3*p*q+c4*q*r+c1*M_x+c2*M_z;
qdot = c7*p*r-c6*(p*p - r*r)+c5*M_y;
rdot = c9*p*q-c3*q*r+c2*M_x+c8*M_z;

Xdot(1:6,1) = [udot vdot wdot pdot qdot rdot];

% phidot = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
% thetadot = q*cos(phi) - r*sin(phi);
% psidot = q*sin(phi)*sec(theta) + r*cos(phi)*sec(theta);

C_euler = [1 sin(phi)*tan(theta) cos(phi)*tan(theta); 
    0 cos(phi) -sin(phi); 
    0 sin(phi)*sec(theta) cos(phi)*sec(theta)];

Xdot(7:9,1) = C_euler * [p;q;r];

% Xdot(1:9,1) = [udot vdot wdot pdot qdot rdot phidot thetadot psidot]';

C_bn = DirectionCosineMatrix(phi,theta,psi);

Xdot(10:12,1) = C_bn' * [u;v;w];

return

function [rho,a] = Atmosphere(alt)

lap_rate= 0.0065;   % (C)deg/meter
t0=15;            % (C)deg
rho0=1.225;         % Kg/m^3
p0=101310;          % Pa


if alt>= 11000,
    temp=-56.4;
    pres=p0*(0.2189*(exp(9.81*(11000-alt)/(287.1*(temp+273)))));
    rho= rho0*(0.2972*(exp(9.81*(11000-alt)/(287.1*(temp+273)))));
else
    temp=(t0-(alt*lap_rate));
    pres=p0*((((temp+273)/(t0+273))^(9.81/(lap_rate*287.1))));
    rho =rho0*((((temp+273)/(t0+273))^((9.81-(lap_rate*287.1))/(lap_rate*287.1))));
end,
a=(1.4*287.1*(temp+273))^0.5;

return