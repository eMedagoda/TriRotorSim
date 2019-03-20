function Xdot = StateRatesRigidBody(X,U)

global Params

g = Params.g;

Ixx = Params.Ixx;
Iyy = Params.Iyy;
Izz = Params.Izz;
Ixz = Params.Ixz;

Ax = Params.Ax;
Ay = Params.Ay;
Az = Params.Az;

m = Params.m;

% calculate body axis rotational inertias
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

u     = X(1);
v     = X(2);
w     = X(3);
p     = X(4);
q     = X(5);
r     = X(6);
phi   = X(7);
theta = X(8);
psi   = X(9);

% % body forces
% rho = Atmosphere(0.0);
% F_x = U(1) - sign(u)*1.28*rho*Ax*0.5*u*u;
% F_y = U(2) - sign(v)*1.28*rho*Ay*0.5*v*v; % added drag to model (Very high lateral drag coef)
% F_z = U(3) - sign(w)*1.28*rho*Az*0.5*w*w;

F_x = U(1);
F_y = U(2);
F_z = U(3);

% body moments
M_x = U(4);
M_y = U(5);
M_z = U(6);

% body accelerations
udot = -g*sin(theta)          + r*v - q*w + F_x/m ;
vdot =  g*sin(phi)*cos(theta) - r*u + p*w + F_y/m ;
wdot =  g*cos(phi)*cos(theta) + q*u - p*v + F_z/m ;

% calculate body axis rotation rates
pdot = c3*p*q + c4*q*r + c1*M_x + c2*M_z;
qdot = c7*p*r - c6*(p*p - r*r)  + c5*M_y;
rdot = c9*p*q - c3*q*r + c2*M_x + c8*M_z;

Xdot(1:6,1) = [udot vdot wdot pdot qdot rdot]';

% % euler angle rates
% phidot = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
% thetadot = q*cos(phi) - r*sin(phi);
% psidot = q*sin(phi)*sec(theta) + r*cos(phi)*sec(theta);

C_euler = [1 sin(phi)*tan(theta) cos(phi)*tan(theta); 
    0 cos(phi) -sin(phi); 
    0 sin(phi)*sec(theta) cos(phi)*sec(theta)];

Xdot(7:9,1) = C_euler * [p;q;r];

% navigation to body
C_bn = DirectionCosineMatrix(phi,theta,psi);

% position rates (navigation frame)
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