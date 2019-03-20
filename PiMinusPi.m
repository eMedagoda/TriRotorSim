function [X] = PiMinusPi(X)

% if vehicle bearing is less than -pi
if X <= -pi

    X = X + 2*pi;

% if vehicle bearing is greater than pi
elseif X >= pi

    X = X - 2*pi;

end

return
