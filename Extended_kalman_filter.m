% Clear workspace and command window
function [xe_total, xt_total] = Extended_kalman_filter(w, v, xtraj, Qw, Rv)
    Q = diag([0.1, 0.1, 0.1, 0.1, 5000, 5000, 100, 100000, 1000000, 100000, 100, 100]);
    R = diag([0.001, 0.01, 0.001, 0.001]);
    step = 0.1; 
    t_max = 0.8; 
    t_current = 0;
    tol = 0.5;  %% very small because angles
    C = eye(12);        % Since we want the output (yn) to be our state: yn = I*xn +[0]un
    D = zeros(12, 4);   % See above
    xe_end = zeros(12, 1); 
    xt_end = zeros(12, 1); 
    t_control = [0]; 
    xe_total = zeros(12, 1); 
    xt_total = zeros(12, 1); 
    for i = 1:size(xtraj, 2)
        x_eq = xtraj(1:end, i) ; 
       
        u_eq = gravity_compensation(x_eq); 
        [Ai, Bi] = jacobian_crane(x_eq, u_eq); 
        K = lqr(Ai, Bi, Q, R); 
        
        kalmf = kalman(ss(Ai, [Bi Bi], C, [D D]), Qw, Rv);
        kalmf = kalmf(13:end,:); %%% disregard the ys because they're just the states*1

        while any(abs(xe_end - x_eq) > tol) 
            t_span = t_current:step:t_current+t_max; 
            [xe_history, xt_history, ~] = Kalman_filter(kalmf, Ai, Bi, C, D, K, t_span, w, v, x_eq, [xe_end; xt_end]);   
            t_control = [t_control t_span]; 
            xe_total = [xe_total xe_history]; 
            xt_total = [xt_total xt_history]; 
            t_current = t_current + t_max; 
            xe_end = xe_history(1:end, end); 
            xt_end = xt_history(1:end, end);

            abs(xe_end - x_eq) 
        end
    
    end
    
    
end 