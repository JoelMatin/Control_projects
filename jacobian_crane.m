function [A, B] = jacobian_crane(x_eq, u_eq)

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
    
    
    MatriX11 = Numerical_matrix(dataM11, x_eq, u_eq); 
    MatriX12 = Numerical_matrix(dataM12, x_eq, u_eq); 
    MatriX13 = Numerical_matrix(dataM13, x_eq, u_eq); 
    MatriX14 = Numerical_matrix(dataM14, x_eq, u_eq); 
    MatriX15 = Numerical_matrix(dataM15, x_eq, u_eq); 
    MatriX16 = Numerical_matrix(dataM16, x_eq, u_eq); 
    MatriX17 = Numerical_matrix(dataM17, x_eq, u_eq); 
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

end