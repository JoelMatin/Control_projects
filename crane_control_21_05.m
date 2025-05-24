close all; clear all; clc; 


%% Crane Control %% 


q_eq = [pi/6; pi/3; -pi/6; 1; 0; 0]; %%%% 1 is the length of cable 

I_tot=250;   l_B=2.5;     m_B=300;     I_B=156.25;     l_J=2;     m_J=250;     I_J=85;     g=9.81;       m=90;

th1 = q_eq(1);
th2 = q_eq(2);
th3 = q_eq(3);
d6  = q_eq(4);
th4 = q_eq(5);
th5 = q_eq(6);

q_eqdot = zeros(6, 1);
dth1 = q_eqdot(1);
dth2 = q_eqdot(2);
dth3 = q_eqdot(3);
dd6  = q_eqdot(4);
dth4 = q_eqdot(5);
dth5  =q_eqdot(6);



%% Linearization %%

u_eq = [                               0;
(g*l_B*cos(th2)*(2*m + m_B + 2*m_J))/2;
     (g*l_J*cos(th3)*(2*m + m_J))/2;
             -g*m*cos(th4)*cos(th5)];           %% gravity compensation term taken from the gravity effect G in the dynamics function but without two last terms because not actuable


x_eq = [zeros(6,1); pi/6; pi/3; -pi/6; 0.5; 0; 0] % [qdot_eq; q_eq] so speeds set to zero because equilibrium point and cable =  0.5m


%[X_eq, u, y, dx] = trim('crane_model_21_05', x_eq, u_eq);   %%%%not useful just to check that we indeed have equilibrium at gravity compensation

[A, B, C, D] = linmod('crane_model_21_05', x_eq, u_eq)    %C and D cannot be computed in this model since our output is the state itself


% C = eye(12);% Since we want the output (yn) to be our state: yn = I*xn +[0]un
% D = zeros(12, 4);   % See above
% No need for this now i've cooked

sys_linearized = ss(A, B, C, D);



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
C_dague = pinv(C)  
B_dague = pinv(B);
N = K*C_dague - B_dague * A * C_dague;  


% Let check system stability:
poles = eig(A)
poles_CL = eig(A - B*K) % OK our lqr stabilized our system :3


Matrix_speed = [eye(6, 12); zeros(6, 12)]
Matrix_pos = diag([zeros(1, floor(12/2)), ones(1, ceil(12/2))])




%% Kalman filter %%
% The aim here is to try to estimate each state depending on the
% observations we do. Since our model doesn't represent the reality.


