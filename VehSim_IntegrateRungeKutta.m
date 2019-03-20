function [X_NEW] = VehSim_IntegrateRungeKutta(X,U,DT)
    
% calculate vehicle state rates
[X_dot1] = VehSim_StateRates(X,U);

% Propagate forward full step
An = X_dot1 * DT;

% calculate vehicle state rates
X_dot2 = VehSim_StateRates(X + 0.5 * An,U);

% Propagate forward 1/2 step
Bn = X_dot2 * DT;  

% calculate vehicle state rates
X_dot3 = VehSim_StateRates(X + 0.5 * Bn,U);

% Update 1/2 step propagation
Cn = X_dot3 * DT;  

% calculate vehicle state rates
X_dot4 = VehSim_StateRates(X + Cn,U);

% Update full step propagation
Dn = X_dot4 * DT;

% Estimate final weighted average state propagation
X_NEW = X + (An + 2*Bn + 2*Cn + Dn)/6;

return