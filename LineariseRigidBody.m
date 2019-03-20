function [A, B] = LineariseRigidBody(X_TRIM,U_TRIM)

Xdot_TRIM = StateRatesRigidBody(X_TRIM,U_TRIM);

num_states = length(X_TRIM);
num_controls = length(U_TRIM);

A = zeros(num_states,num_states);
B = zeros(num_states,num_controls);

dX = 1e-7;

for i = 1:num_states

    X_PERT = X_TRIM;

    X_PERT(i) = X_PERT(i) + dX;

    Xdot_PERT_X = StateRatesRigidBody(X_PERT,U_TRIM);

    A(:,i) = (Xdot_PERT_X - Xdot_TRIM)./dX;

end

for i = 1:num_controls

    U_PERT = U_TRIM;

    U_PERT(i) = U_PERT(i) + dX;

    Xdot_PERT_U = StateRatesRigidBody(X_TRIM,U_PERT);

    B(:,i) = (Xdot_PERT_U - Xdot_TRIM)./dX;

end

return