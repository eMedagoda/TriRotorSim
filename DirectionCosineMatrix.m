function [C_bn] = DirectionCosineMatrix(phi,theta,psi)

% rotational transformation about z-axis
C3 = L3(psi);

% rotational transformation about y-axis
C2 = L2(theta);

% rotational transformation about x-axis
C1 = L1(phi);

% transformation matrix from navigation frame to body frame
C_bn = [C1*C2*C3];

return