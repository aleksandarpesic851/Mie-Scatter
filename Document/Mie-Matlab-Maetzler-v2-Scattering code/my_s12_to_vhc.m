 v = 30;
 lamda = 0.842;
 ns = 1.33;
 n0 = 1.335;
 alpha = 0.001942;
 hc = 31;
 m_mass = 66500;
 eps = 0.252;
 m1 = n0 + alpha * hc;
 m2 = log(10) * lamda * eps * hc / (pi * m_mass);
 m = (m1 - 1i*m2) / ns;
 
 r = ((3 * v) / (4 * pi))^(1/3);
 k = 2 * pi * ns / lamda;
 x = k * r;
 res = [];
 
 for theta = 0:0.1:21
     u = pi * theta / 90;
     cosU = cos(u);
     
     s12 = Mie_S12(m, x, cosU);
     s1 = s12(1);
     s2 = s12(2);
     
     resS = (abs(s1) ^ 2 + abs(s2) ^ 2) * lamda ^ 2 * sin(u) / ( 8 * pi ^ 2 * ns ^ 2);
     res = [res, resS];
 end
 theta = 0:0.1:21;
 plot(theta, res);