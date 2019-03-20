function [L1] = L1(phi)

L1 = [1        0         0  ;
    0  cos(phi)  sin(phi) ;
    0  -sin(phi) cos(phi)];

return