function [Synthetic] = SyntheticWaypointDynamics(X,V,NAV,Waypoint,Synthetic,DT)

% calculate current virtual range to virtual target
Synthetic.Range = sqrt((Synthetic.Waypoint(1,NAV) - X(1))^2 + ...
    (Synthetic.Waypoint(2,NAV) - X(2))^2);

% determine virtual range
VirtualRange = Synthetic.TimeHorizon * V;

% determine the desired target velocity
V_w = V * (VirtualRange/Synthetic.Range);

% calculate current bearing to waypoints
Synthetic.PsiC = atan2((Synthetic.Waypoint(2,NAV) - X(2)), ...
    (Synthetic.Waypoint(1,NAV) - X(1)));

% maintain heading within bounds
Synthetic.PsiC = PiMinusPi(Synthetic.PsiC);

% determine new target position to maintain virtual
% range
Synthetic.Waypoint(1:2,NAV) = Synthetic.Waypoint(1:2,NAV) ...
    + V_w * [sin(Waypoint(3,NAV)) cos(Waypoint(3,NAV))]' * DT;

return