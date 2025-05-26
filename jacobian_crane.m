
I_tot=250;   l_B=2.5;     m_B=300;     I_B=156.25;     l_J=2;     m_J=250;     I_J=85;     g=9.81;       m=90;

u_eq = [                               0;
(g*l_B*cos(th2)*(2*m + m_B + 2*m_J))/2;
     (g*l_J*cos(th3)*(2*m + m_J))/2;
             -g*m*cos(th4)*cos(th5)];           %% gravity compensation term taken from the gravity effect G in the dynamics function but without two last terms because not actuable


x_eq = [zeros(6,1); pi/6; pi/3; -pi/6; 0.5; 0; 0] % [qdot_eq; q_eq] so speeds set to zero because equilibrium point and cable =  0.5m


dataM11 = readcell('M11.csv'); 
dataM12 = readcell('M12.csv'); 
dataM13 = readcell('M13.csv'); 
dataM14 = readcell('M14.csv'); 
dataM15 = readcell('M15.csv'); 
dataM16 = readcell('M16.csv'); 
dataM17 = readcell('M17.csv'); 
dataM18 = readcell('M18.csv'); 
dataM19 = readcell('M19.csv'); 
dataM110 = readcell('M110.csv'); 
dataM111 = readcell('M111csv.csv'); 
dataM112 = readcell('M112.csv'); 
dataB11 = readcell('B11.csv'); 
dataB12 = readcell('B12.csv'); 
dataB13 = readcell('B13.csv'); 
dataB14 = readcell('B14.csv'); 


MatriX11 = zeros(6, 1); %Numerical_matrix(dataM11, x_eq, u_eq); 
MatriX12 = Numerical_matrix(dataM12, x_eq, u_eq); 
MatriX13 = Numerical_matrix(dataM13, x_eq, u_eq); 
MatriX14 = Numerical_matrix(dataM14, x_eq, u_eq); 
MatriX15 = Numerical_matrix(dataM15, x_eq, u_eq); 
MatriX16 = Numerical_matrix(dataM16, x_eq, u_eq); 
MatriX17 = zeros(6, 1); %Numerical_matrix(dataM17, x_eq, u_eq); 
MatriX18 = Numerical_matrix(dataM18, x_eq, u_eq); 
MatriX19 = Numerical_matrix(dataM19, x_eq, u_eq); 
MatriX110 = Numerical_matrix(dataM110, x_eq, u_eq); 
MatriX111 = Numerical_matrix(dataM111, x_eq, u_eq); 
MatriX112 = Numerical_matrix(dataM112, x_eq, u_eq); 

B11 = Numerical_matrix(dataB11, x_eq, u_eq); 
B12 = Numerical_matrix(dataB12, x_eq, u_eq); 
B13 = Numerical_matrix(dataB13, x_eq, u_eq); 
B14 = Numerical_matrix(dataB14, x_eq, u_eq); 

A = [MatriX11, MatriX12, MatriX13, MatriX14, MatriX15, MatriX16,  MatriX17, MatriX18, MatriX19, MatriX110, MatriX111, MatriX112; 
     eye(6, 6), zeros(6, 6)]; 
 
B = [B11, B12, B13, B14; 
    zeros(6, 4) ]; 