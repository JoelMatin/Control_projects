% Clear workspace and command window
function [xe_history, xt_history, u_history] = Kalman_filter(kalmf, A, B, C, D, K, tspan, w, v, xd, z0)
    
    n_states = size(A, 1);
    n_inputs = size(B, 2);

    % Extract Kalman filter matrices
    A_kf = kalmf.A;
    B_kf = kalmf.B;
    
    % Define ODE function for plant and Kalman filter
    function dzdt = odefun(t, z, A, B, C, D, A_kf, B_kf, K, w, v, xd, tspan, gravity_compensation)
        % Combined state vector: z = [xt; xe]
        xt = z(1:12); % True states
        xe = z(13:24); % Estimated states
        
        % Interpolate noise at time t
        idx = max(1, min(floor(t / (tspan(2) - tspan(1))) + 1, length(tspan)));
        w_t = w(:, idx);
        v_t = v(:, idx);
        
        % LQR controller
        xe_tilde = xe - xd;
        u_tilde = -K * xe_tilde;
        
        % Gravity compensation
        ud = gravity_compensation(xe);
        
        % Total control input:
        u = u_tilde + ud;
        
        % Measurements out of the plant
        x = C * xt + [D D] * [u; w_t] + v_t; %% noisy states
        
        % Disturbed plant dynamics
        xt_dot = A * xt + [B B] * [u; w_t];
        
        % Kalman filter dynamics
        kf_input = [u; x];
        xe_dot = A_kf * xe + B_kf * kf_input;
        
        % Combined dynamics
        dzdt = [xt_dot; xe_dot];
    end
    
    % Initial conditions: [xt; xe]
    % z0 = zeros(2 * n_states, 1); % Initial true and estimated states

    % Solve ODE using ode15
    [t, z] = ode15s(@(t,z) odefun(t, z, A, B, C, D, A_kf, B_kf, K, w, v, xd, tspan, @gravity_compensation), tspan, z0);
    
    % Extract results
    xt_history = z(:,1:12)'; % True states
    xe_history = z(:,13:24)'; % Estimated states
    
    % Compute control inputs
    u_history = zeros(n_inputs, length(t));
    for k = 1:length(t)
        xe_tilde = xe_history(:,k) - xd;
        u_tilde = -K * xe_tilde;
        ud = gravity_compensation(xe_history(:,k));
        u_history(:,k) = u_tilde + ud;
    end
    
end 