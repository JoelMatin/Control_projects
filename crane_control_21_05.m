close all; clear all; clc; 

%% Crane Control %%% 


I_tot=250;   l_B=2.5;     m_B=300;     I_B=156.25;     l_J=2;     m_J=250;     I_J=85;     g=9.81;       m=90;

x_eq = [zeros(6,1); pi/6; pi/3; -pi/6; 0.5; 0; 0] ; % [qdot_eq; q_eq] so speeds set to zero because equilibrium point and cable =  0.5m

th1 = x_eq(7);
th2 = x_eq(8);
th3 = x_eq(9);
d6  = x_eq(10);
th4 = x_eq(11);
th5 = x_eq(12);

u_eq = [                               0;
(g*l_B*cos(th2)*(2*m + m_B + 2*m_J))/2;
     (g*l_J*cos(th3)*(2*m + m_J))/2;
             -g*m*cos(th4)*cos(th5)];           %% gravity compensation term taken from the gravity effect G in the dynamics function but without two last terms because not actuable


%% Linearization using linmod %%

%[X_eq, u, y, dx] = trim('crane_model_21_05', x_eq, u_eq);   %%%%not useful just to check that we indeed have equilibrium at gravity compensation

[Alin, Blin, C, D] = linmod('crane_model_21_05', x_eq, u_eq)    %C and D cannot be computed in this model since our output is the state itself

sys_linearized = ss(Alin, Blin, C, D);

%% Linearization using Jacobian analxtic derivation

% load("matrices_jacobian_better.mat")
[A, B] = jacobian_crane(x_eq, u_eq); 
C = eye(12);        % Since we want the output (yn) to be our state: yn = I*xn +[0]un
D = zeros(12, 4);   % See above
sys_linearized_jac = ss(A, B, C, D); 

%Indeed both ways give the same matrices
A-Alin
B-Blin

%% Reference 
% Moving a little bit compared to the eq values since we are linearizing
% around this eq point

x_d = [0;
       0;
       0;
       0;
       0;
       0;
      pi/10;
      pi/10;
      pi/10;
      0.01;
      0;
      0] + x_eq

%% LQR Control %% 
% Lets now do a full state feedback with feedforward term for ref tracking
% i.e u = -K*x + N*r. Please see CSD course for ref.
% The general idea: x_dot = 0 at steady state and y = C*x = r 

x0 = x_eq;
diagQ = ones(12,1);         %% Multiplies the states so 12x12 
diagR = ones(4,1);          %% Multiplies input so 4x4
Q = diag(diagQ);
R = diag(diagR);

% Tuning Q and R matrixes to assess good performances: this tuning shows we
% really want to penalize the behavior and not the input

Q = diag([0.1, 0.1, 0.1, 0.1, 5000, 5000, 100, 100000, 1000000, 100000, 100, 100]);
R = diag([0.001, 0.01, 0.001, 0.001]);

K = lqr(A, B, Q, R);

% Lets now compute N
C_dague = pinv(C)  ; 
B_dague = pinv(B);
N = K*C_dague - B_dague * A * C_dague;  


% Let check system stabili ty:
poles = eig(A); 
poles_CL = eig(A - B*K);  % OK our lqr stabilized our system :3

%%% Matrices to separate the speed states from the position states in the
%%% simulink scopes
Matrix_speed = [eye(6, 12); zeros(6, 12)]; 
Matrix_pos = diag([zeros(1, floor(12/2)), ones(1, ceil(12/2))]); 

%%% Results: LQR indeed controls the crane both in linearized (with and
%%% without jacobian ) and in real model but we have: steady state error,
%%% oscillations for both posiitons and speed and a long rising time for
%%% the q1



%% Kalman filter %%
% The aim here is to try to estimate each state depending on the
% observations we do. Since our model doesn't represent the reality.
clc; 

% %%% Covariance of input process disturbance and output sensor noise
Q = 1 * eye(4); 
R = 1 * eye(12);  

% plant_disturbed = ss(A, [B B], C, [D D]); 

[kalmf, L, P, Mx, Z] = kalman(ss(A, [B B], C, [D D]), Q, R);
kalmf = kalmf(13:end,:); %%% disregard the ys because they're just the states*1

% plant_disturbed.InputName = {'u(1)', 'u(2)', 'u(3)', 'u(4)', 'w1', 'w2', 'w3', 'w4'};
% plant_disturbed.OutputName = {'xt(1)', 'xt(2)', 'xt(3)', 'xt(4)', 'xt(5)', 'xt(6)', 'xt(7)', 'xt(8)', 'xt(9)', 'xt(10)', 'xt(11)', 'xt(12)'};
% vIn = sumblk('x = xt + v', 12);         %% define output corrupted by noise 
% 

kalmf.InputName = {'u(1)', 'u(2)', 'u(3)', 'u(4)', 'x(1)', 'x(2)', 'x(3)', 'x(4)', 'x(5)', 'x(6)', 'x(7)', 'x(8)', 'x(9)', 'x(10)', 'x(11)', 'x(12)'};
kalmf.OutputName = {'xe(1)', 'xe(2)', 'xe(3)', 'xe(4)', 'xe(5)', 'xe(6)', 'xe(7)', 'xe(8)', 'xe(9)', 'xe(10)', 'xe(11)', 'xe(12)'};    %% estimated states 

% model = 'Kalman_filter_crane'; %to have the gravity compensation% Add MATLAB Function block for gravity compensation (if not already added)
% gravity_comp_block = [model '/gravity_comp']; % Reference the block by its path
% 
% xdIn = sumblk('xe_tilde = xe - xd', 12);  %% 
% % Define LQR controller: u = -K*xe + N*xd; actually this was before but
% % now i'm doing u = -K*utilde+ud 
% lqr_controller = ss(-K); % Static gain for u_tilde = -K*xe_tilde
% lqr_controller.InputName = {'xe_tilde(1)', 'xe_tilde(2)', 'xe_tilde(3)', 'xe_tilde(4)', 'xe_tilde(5)', 'xe_tilde(6)', ...
%                             'xe_tilde(7)', 'xe_tilde(8)', 'xe_tilde(9)', 'xe_tilde(10)', 'xe_tilde(11)', 'xe_tilde(12)'}; 
% lqr_controller.OutputName = {'u_tilde(1)', 'u_tilde(2)', 'u_tilde(3)', 'u_tilde(4)'};
% 
% udIn = sumblk('u = u_tilde + ud', 4);
% 
% % Connect the system
% plant_filtered = connect(plant_disturbed, vIn, kalmf, xdIn, lqr_controller, gravity_comp_block, udIn, ...
%     {'w1', 'w2', 'w3', 'w4', 'v(1)', 'v(2)', 'v(3)', 'v(4)', 'v(5)', 'v(6)', ...
%      'v(7)', 'v(8)', 'v(9)', 'v(10)', 'v(11)', 'v(12)', ...
%      'xd(1)', 'xd(2)', 'xd(3)', 'xd(4)', 'xd(5)', 'xd(6)', ...
%      'xd(7)', 'xd(8)', 'xd(9)', 'xd(10)', 'xd(11)', 'xd(12)'}, ...
%     {'xt(1)', 'xt(2)', 'xt(3)', 'xt(4)', 'xt(5)', 'xt(6)', ...
%      'xt(7)', 'xt(8)', 'xt(9)', 'xt(10)', 'xt(11)', 'xt(12)', ...
%      'xe(1)', 'xe(2)', 'xe(3)', 'xe(4)', 'xe(5)', 'xe(6)', ...
%      'xe(7)', 'xe(8)', 'xe(9)', 'xe(10)', 'xe(11)', 'xe(12)'});

%% Input Signals generation

t_vec = (0:0.1:50)';
n = length(t_vec);
x_eq_traj = zeros(12, n);

q_eq1 = zeros(6, 1); 
q_eq2 = x_eq(7:end); %[pi/6; pi/3; -pi/6; 0.5; 0; 0];
q_eq3 = [pi/3; pi/6; pi/4; 1; 0; 0]; 
% q_eq4 = [pi/6; pi/3; pi/3; 0.5; 0; 0]; 
q_eqs = [q_eq1, q_eq2]; 

% Time points for equilibrium configurations
t_eqs = zeros(1, size(q_eqs, 2)); 
for i = 2:length(t_eqs)
    t_eqs(i) = 0 + (i)*max(t_vec)/length(t_eqs); 
end

% x_eq_traj(7, :) = interp1(t_eqs, q_eqs(1,:), t_vec, 'linear');
% x_eq_traj(8, :) = interp1(t_eqs, q_eqs(2,:), t_vec, 'linear');
% x_eq_traj(9, :) = interp1(t_eqs, q_eqs(3,:), t_vec, 'linear');
% x_eq_traj(10, :) = interp1(t_eqs, q_eqs(4,:), t_vec, 'linear');
% x_eq_traj(11, :) = interp1(t_eqs, q_eqs(5,:), t_vec, 'linear');
% x_eq_traj(12, :) = interp1(t_eqs, q_eqs(6,:), t_vec, 'linear');

for i = 7:size(x_eq_traj, 1)
    x_eq_traj(i, :) = q_eq2(i-6).*ones(1, 501);
end  


% Generate process and sensor noise (piecewise constant for ode15)
rng(0); % For reproducibility
w = sqrt(Q) * randn(n_noise, n_steps); % Process noise
v = sqrt(R) * randn(n_meas, n_steps); % Sensor nois
xd = x_eq_traj; 


% % Plotting
% figure(2)
% 
% % Subplot for x_eq (configuration states)
% % subplot(4, 1, 1);
% plot(t_vec, x_eq_traj(7, :), 'b-', 'LineWidth', 1.5, 'DisplayName', '\theta_1');
% hold on;
% plot(t_vec, x_eq_traj(8, :), 'r-', 'LineWidth', 1.5, 'DisplayName', '\theta_2');
% plot(t_vec, x_eq_traj(9, :), 'g-', 'LineWidth', 1.5, 'DisplayName', '\theta_3');
% plot(t_vec, x_eq_traj(10, :), 'k-', 'LineWidth', 1.5, 'DisplayName', 'cable');hold on; 
% plot(t_vec, x_eq_traj(11, :), 'm-', 'LineWidth', 1.5, 'DisplayName', '\theta_4');
% plot(t_vec, x_eq_traj(12, :), 'c-', 'LineWidth', 1.5, 'DisplayName', '\theta_5');
% hold off; 
% title('Equilibrium State Trajectory');
% xlabel('Time (s)');
% ylabel('State Value');
% legend('show', 'Location', 'best');
% grid on;


%plot noise
% 
% % Subplot 3: Process noise w
% subplot(4, 1, 3);
% plot(t_vec, w(1, :), 'b-', 'LineWidth', 1.5, 'DisplayName', 'w_1');
% hold on;
% plot(t_vec, w(2, :), 'r-', 'LineWidth', 1.5, 'DisplayName', 'w_2');
% plot(t_vec, w(3, :), 'g-', 'LineWidth', 1.5, 'DisplayName', 'w_3');
% plot(t_vec, w(4, :), 'k-', 'LineWidth', 1.5, 'DisplayName', 'w_4');
% hold off;
% title('Process Noise');
% xlabel('Time (s)');
% ylabel('Noise Amplitude');
% legend('show', 'Location', 'best');
% grid on;
% 
% % Subplot 4: Measurement noise v (plotting first 4 components to avoid clutter)
% subplot(4, 1, 4);
% plot(t_vec, v(1, :), 'b-', 'LineWidth', 1.5, 'DisplayName', 'v_1');
% hold on;
% plot(t_vec, v(2, :), 'r-', 'LineWidth', 1.5, 'DisplayName', 'v_2');
% plot(t_vec, v(3, :), 'g-', 'LineWidth', 1.5, 'DisplayName', 'v_3');
% plot(t_vec, v(4, :), 'k-', 'LineWidth', 1.5, 'DisplayName', 'v_4');
% hold off;
% title('Measurement Noise (First 4 Components)');
% xlabel('Time (s)');
% ylabel('Noise Amplitude');
% legend('show', 'Location', 'best');
% grid on;


%% Simulation with ode15s 

[xe_history, xt_history, u_history] = Kalman_filter(kalmf, A, B, C, D, K, t_vec, w, v, x_eq); 


% Plot results ith state
figure;
subplot(3,1,1);
i = 7
plot(t_vec, xt_history(i,:), 'b-', 'DisplayName', 'True State (xt(i))');
hold on;
plot(t_vec, xd(i,:), 'b-', 'DisplayName', 'Reference State (xt(i))');
hold on;
plot(t_vec, xe_history(i,:), 'r--', 'DisplayName', 'Estimated State (xe(i))');
xlabel('Time (s)');
ylabel('State 1');
legend;
title('True vs Estimated States');

subplot(3,1,2);
plot(t_vec, xt_history(i,:) - xe_history(i,:), 'k-', 'DisplayName', 'Estimation Error');
xlabel('Time (s)');
ylabel('Error (xt(1) - xe(1))');
legend;

subplot(3,1,3);
plot(t_vec, u_history(1,:), 'g-', 'DisplayName', 'Control Input (u(i))');
xlabel('Time (s)');
ylabel('Control Input 1');
legend;

%% Simulation with lsim and simulink
% Combine inputs: [w; v; xd]
% U = [w; v; x_eq_traj]; % 28xn input vector: [4 process noise + 12 measurement noise + 12 reference]

% Initial condition
x0 = zeros(12, 1); % Initial state for plant
x0_kalman = zeros(12, 1); % Initial state for Kalman filter
X0 = [x0; x0_kalman]; % Combined initial condition

% Simulate the system
% [Y, T, X] = lsim(plant_filtered, U, t_vec, X0);

v = [t_vec, v']; 
w = [t_vec, w']; 
xd = [t_vec, xd']; 

sim('Kalman_filter_crane'); 

% Extract true and estimated states
xe = xestimated.Data; % Estimated states (xe(1) to xe(12))
xt = xt.Data; % Estimated states (xe(1) to xe(12))
xd = xdout.Data'; % Estimated states (xe(1) to xe(12))
tout = xestimated.Time; 

% Plotting Results
figure(6);
i = 7;
plot(tout, 180*xt(i, :)/pi, 'b-', 'LineWidth', 1.5, 'DisplayName', 'True');
hold on;
plot(tout, 180*xe(i, :)/pi, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Estimated');
plot(tout, 180*xd(i, :)/pi, 'g:', 'LineWidth', 1, 'DisplayName', 'Reference');
hold off;
title(['State x(', num2str(i), ')']);
xlabel('Time (s)');
ylabel(['x_', num2str(i), '(degrees)']);
legend('show');
grid on;


