function [r0, v0] = Kep2Car(a, e, i, Raan, omega, theta, mu)

p = a * (1 - e^2);
v = sqrt(mu/p);
r = p/(1+ e*cos(theta));

% definisco nel sistema perifocale i valori di raggio e velcità

r_pqw = [r*cos(theta) r*sin(theta) 0]';
v_pqw =  v * [-sin(theta) e+cos(theta) 0]';

% la direzione verticale mi è completamente data dall'inclinazione del
% piano orbitale


R = [ cos(Raan)*cos(omega) - sin(Raan)*sin(omega)*cos(i) -cos(Raan)*sin(omega) - sin(Raan)*cos(omega)*cos(i) sin(Raan)*sin(i);
      sin(Raan)*cos(omega) + cos(Raan)*sin(omega)*cos(i) -sin(Raan)*sin(omega) + cos(Raan)*cos(omega)*cos(i) -cos(Raan)*sin(i);
      sin(omega)*sin(i), cos(omega)*sin(i), cos(i)];

% Coordinate nel sistema inerziale

r0 = R * r_pqw;
v0 = R * v_pqw;

end