 hc = 20; 
 lamda = 0.6328;
 ns = 1.3339;
 n0 = 1.345;
 alpha = 0.001942;
 m_mass = 66500;
 eps = 0.1616;
 m1 = n0 + alpha * hc;
 m2 = log(10) * lamda * eps * hc / (pi * m_mass);
 m = (m1 - 1i*m2) / ns;
 k = 2 * pi * ns / lamda;
 
 
 res = [];
 
 theta = 5;
 d_theta = 15;
 theta_step = 0.01;
 
 for v = 30:180
     
     r = ((3 * v) / (4 * pi))^(1/3);
     x = k * r;
     integ_s = 0;
     for beta = theta:theta_step:theta+d_theta
         u = pi * beta / 180;
         cosU = cos(u);

         s12 = Mie_S12(m, x, cosU);
         s1 = s12(1);
         s2 = s12(2);

         integ_s = integ_s + (abs(s1) ^ 2 + abs(s2) ^ 2) * lamda ^ 2 * sin(u) / ( 8 * pi ^ 2 * ns ^ 2);
     end
     integ_s = integ_s * theta_step * pi / 180;
     
     res = [res, integ_s];
 end
 
 v = 30:180;
 plot(v, res);