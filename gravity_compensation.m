function u_comp = gravity_compensation(x)

I_tot=250;   l_B=2.5;     m_B=300;     I_B=156.25;     l_J=2;     m_J=250;     I_J=85;     g=9.81;       m=90;
  
q = x(7:12);

th1 = q(1);
th2 = q(2);
th3 = q(3);
d6  = q(4);
th4 = q(5);
th5 = q(6);
    
u_comp = [                               0;
(g*l_B*cos(th2)*(2*m + m_B + 2*m_J))/2;
     (g*l_J*cos(th3)*(2*m + m_J))/2;
             -g*m*cos(th4)*cos(th5)];

end