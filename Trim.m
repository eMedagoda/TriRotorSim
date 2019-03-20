function [X0,U0,n] = Trim(VTrim,AltTrim,ThetaTrim)

% initial trim vector
XTrim = [1 2 3]'; % Lift, Thrust, Pitch

% convergence tolerance
Tol = 1e-10;

% perturbation
dXTrim = 1e-8;

% intial state vector
X = zeros(12,1);

% initial control vector
U = zeros(5,1);

% initialise jacobian
J = zeros(3,3);

% intial error flag
Err = 1;

% initial counter
n = 1;

while Err > Tol
    
    X(1) = VTrim*cos(ThetaTrim);
    X(3) = VTrim*sin(ThetaTrim);
    X(8) = ThetaTrim;
    
    U(1) = XTrim(1);
    U(2) = XTrim(2);
    U(4) = XTrim(3);
    
    Xdot = TriSim_StateRates(X,U);
    
    XTrimDot = [Xdot(1) Xdot(3) Xdot(5)]';
           
    % perturb controls
    Pert = XTrim(1) + dXTrim;
    
    U(1) = Pert;
    U(2) = XTrim(2);
    U(4) = XTrim(3);
    
    Xdot = TriSim_StateRates(X,U);
    
    XTrimDotPert = [Xdot(1) Xdot(3) Xdot(5)]';
    
    J(:,1) = (XTrimDotPert - XTrimDot)/dXTrim;
    
    % perturb controls
    Pert = XTrim(2) + dXTrim;
    
    U(1) = XTrim(1);
    U(2) = Pert;
    U(4) = XTrim(3);
    
    Xdot = TriSim_StateRates(X,U);
    
    XTrimDotPert = [Xdot(1) Xdot(3) Xdot(5)]';
    
    J(:,2) = (XTrimDotPert - XTrimDot)/dXTrim;
    
    % perturb controls
    Pert = XTrim(3) + dXTrim;
    
    U(1) = XTrim(1);
    U(2) = XTrim(2);
    U(4) = Pert;
    
    Xdot = TriSim_StateRates(X,U);
    
    XTrimDotPert = [Xdot(1) Xdot(3) Xdot(5)]';
    
    J(:,3) = (XTrimDotPert - XTrimDot)/dXTrim;
    
    % check for convergence
    XTrim_new = XTrim - inv(J) * XTrimDot;
    
    Err = abs(XTrim_new - XTrim);
    
    XTrim = XTrim_new;
    
    n = n + 1;
    
end

X0 = zeros(12,1);
U0 = zeros(5,1);

X0(1) = VTrim*cos(ThetaTrim);
X0(3) = VTrim*sin(ThetaTrim);
X0(8) = ThetaTrim;
X0(12) = AltTrim;

U0(1) = XTrim_new(1);
U0(2) = XTrim_new(2);
U0(4) = XTrim_new(3);

return
