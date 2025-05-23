close all; clear all; clc; 


%% LQR with partial feedback linearization %% 


q_eq = [pi/6; pi/3; pi/6; 1; 0; 0]; 

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



%% Pt of equilibrum  and initial state%%

u_eq = [                               0;
(g*l_B*cos(th2)*(2*m + m_B + 2*m_J))/2;
     (g*l_J*cos(th3)*(2*m + m_J))/2;
             -g*m*cos(th4)*cos(th5)];


x_eq = [zeros(6,1); pi/6; pi/3; -pi/6; 0.5; 0; 0] % [qdot_eq; q_eq]

x0 = x_eq;

%% Reference 


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
      0] + x_eq;

xa_d = [x_d(1:4);x_d(7:10)];





%% Outer loop 
%  We have xa_dot = Axa + Bnu  (LINEAR system thanks to the inner loop), let derive A and B to find
%  our LQR K gains. xa = [qa_dot; qa]. Since we have qa_ddot = nu,
%  it follows that xa_dot = [qa_ddot; qa_dot] = [nu, qa_dot] 
%  The A and B matrices derivation is straight forward A = [0 0; I_4x4 0]
%  B = [I_4x4; 0_4X4]

Aa = [zeros(4,4), zeros(4,4);
     eye(4,4)  , zeros(4,4)];

Ba = [eye(4,4);
     zeros(4,4)];

% LQR control law nu = -Ka * xa_tilde, Q is an 8x8 and R a 4x4

Qa = diag([0.1, 0.1, 0.1, 0.1, 100, 100000, 1000000, 100000]);
Ra = diag([0.001, 0.01, 0.001, 0.001]);

Ka = lqr(Aa, Ba, Qa, Ra);




%% Random
% 
% C_dague = pinv(C);  
% B_dague = pinv(B);
% N = K*C_dague - B_dague * A * C_dague;  


% Let check system stability:
poles = eig(Aa)
poles_CL = eig(Aa - Ba*Ka) % OK our lqr stabilized our system :3




