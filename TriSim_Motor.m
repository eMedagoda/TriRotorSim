function [H,n] = TriSim_Motor(U0)

% initial trim vector
XTrim = [13.4 13.4 2.5 0.0 0.0]'; % T_R, T_L, T_T, mu, d_mu (needs to be close to level trim values)

% convergence tolerance
Tol = 1e-6;

% perturbation
dXTrim = 1e-6;

% initial control vector
H = zeros(5,1);

% initialise jacobian
J = zeros(5,5);

% intial error flag
Err = 1;

% initial counter
n = 1;

while Err > Tol
    
    U = ControlMap(XTrim);
    
    dU = U0 - U;
    
    for i = 1:5
        
        % perturb controls
        XTrim_Pert = XTrim;
        
        XTrim_Pert(i) = XTrim(i) + dXTrim;
        
        U_Pert = ControlMap(XTrim_Pert);
        
        dU_Pert = U0 - U_Pert;
        
        J(:,i) = (dU_Pert - dU)/dXTrim;
        
    end        
    
    % check for convergence
    XTrim_new = XTrim - inv(J) * dU;
    
    Err = abs(XTrim_new - XTrim);
    
    XTrim = XTrim_new;
    
    n = n + 1;
    
end

H = XTrim_new;

return

function [U] = ControlMap(H)

global Params

U = zeros(5,1);

U(1) = -H(1)*cos(H(4) + H(5)) - H(2)*cos(H(4) - H(5)) - H(3);

U(2) = H(1)*sin(H(4) + H(5)) + H(2)*sin(H(4) - H(5));

U(3) = -H(1)*Params.L2*cos(H(4) + H(5)) + H(2)*Params.L4*cos(H(4) - H(5));

U(4) = H(1)*Params.L1*cos(H(4) + H(5)) + H(2)*Params.L3*cos(H(4) - H(5)) - H(3)*Params.L5;

U(5) = -H(1)*Params.L2*sin(H(4) + H(5)) + H(2)*Params.L4*sin(H(4) - H(5));

return
