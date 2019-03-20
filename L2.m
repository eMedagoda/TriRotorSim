function [L2] = L2(theta)

L2 = [cos(theta) 0  -sin(theta) ;
    0     1       0      ;
    sin(theta) 0   cos(theta)];

return