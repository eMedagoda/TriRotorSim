function [L3] = L3(psi)

L3 = [ cos(psi)  sin(psi) 0 ;
    -sin(psi)  cos(psi) 0 ;
    0         0     1];

return