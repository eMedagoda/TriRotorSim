function [Xdot] = VehSim_StateRates(X,U)

Xdot = zeros(3,1);

% rear wheel drive, front wheel steer
Xdot(1) = U(1)*cos(X(3));
Xdot(2) = U(1)*sin(X(3));
Xdot(3) = (U(1)/3)*tan(U(2));

return