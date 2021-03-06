function [V_C, psi_C, R] = TriSim_Guidance(P_V,P_T)

x_V = P_V(1);
y_V = P_V(2);

x_T = P_T(1);
y_T = P_T(2);

V_max = 2;
k_v = 1.5;

r_x = x_T - x_V;
r_y = y_T - y_V;

psi_C = atan2(r_y,r_x);

psi_C = PiMinusPi(psi_C);

R = sqrt(r_x*r_x + r_y*r_y);

V_C = (2*V_max/pi)*atan(k_v * (R-1));

return