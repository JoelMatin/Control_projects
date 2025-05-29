function Matrix = Numerical_matrix(data, x_eq, u_eq)   

    I_tot = 250; l_B = 2.5; m_B = 300; I_B = 156.25; l_J = 2; m_J = 250; I_J = 85; g = 9.81; m = 90;

    % Step 2: Define symbolic variable values
    dth1 = x_eq(1);
    dth2 = x_eq(2);
    dth3 = x_eq(3);
    dd6  = x_eq(4);
    dth4 = x_eq(5); 
    dth5 = x_eq(6); 
    th1 = x_eq(7);
    th2 = x_eq(8);
    th3 = x_eq(9);
    d6  = x_eq(10);
    th4 = x_eq(11); 
    th5 = x_eq(12); 

    u1 = u_eq(1);
    u2 = u_eq(2);
    u3 = u_eq(3);
    u4  = u_eq(4);

    % Step 3: Convert to MATLAB syntax and evaluate
    Matrix = zeros(6,1);  % Preallocate 6x1 numeric matrix

    for i = 1:6
        expr = data{i};                            % e.g., 'Sin[th1_eq]'
        if isnumeric(expr)
             Matrix(i) = expr; 
        else
            expr = strrep(expr, '[', '(');            % Convert brackets
            expr = strrep(expr, ']', ')');
            expr = strrep(expr, 'Sin', 'sin');        % Convert function names
            expr = strrep(expr, 'Cos', 'cos');
            expr = strrep(expr, 'Tan', 'tan');
            expr = strrep(expr, 'Pi', 'pi');          % Just in case
            Matrix(i) = eval(expr);          % Evaluate expression numerically
        end
    end

end
