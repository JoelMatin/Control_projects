function x_out = my_state_transition_fcn(x,u)

% Input : states at time step k
% Output : states at time step k+1
Ts = 1/1000;
q_ddot = crane_model_final(u, x);       % Compute q_ddot (x_dot = f(x,u))

q_dot = x(1:6) + q_ddot * Ts;           % Integrate using euler integration
q = x(7:12) + q_dot * Ts;               % Same but for pos 

x_out = [q_dot; q];                     % State vector at time k+1
end