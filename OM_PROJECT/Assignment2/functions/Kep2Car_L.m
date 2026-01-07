function [r0, v0] = Kep2Car_L(a, e, i, Raan, omega, theta, mu)

%KEP2CAR Conversion from Keplerian elements to Cartesian state vector
%
%   [r0, v0] = Kep2Car(a, e, i, Raan, omega, theta, mu)
%
% PROTOTYPE
%   [r0, v0] = Kep2Car(a, e, i, Raan, omega, theta, mu)
%
% DESCRIPTION
%   This function converts a set of classical Keplerian orbital elements,
%   provided as individual scalar inputs, into the corresponding Cartesian
%   position and velocity vectors expressed in an inertial reference frame
%   (ECI).
%   The conversion is performed by first computing the position and velocity
%   in the perifocal reference frame (PQW) and then rotating them into the
%   inertial frame through a standard 3-1-3 Euler rotation sequence.
%
% INPUT
%   a      [1x1] Semi-major axis                      [km]
%   e      [1x1] Eccentricity                         [-]
%   i      [1x1] Inclination                          [rad]
%   Raan   [1x1] Right ascension of the ascending node[rad]
%   omega  [1x1] Argument of pericenter               [rad]
%   theta  [1x1] True anomaly                         [rad]
%   mu     [1x1] Gravitational parameter of the
%                central body                         [km^3/s^2]
%
% OUTPUT
%   r0 [3x1] Position vector in inertial frame (ECI)  [km]
%   v0 [3x1] Velocity vector in inertial frame (ECI)  [km/s]
%
% ASSUMPTIONS
%   - Two-body Keplerian motion
%   - Elliptic orbit (0 <= e < 1)
%   - Inertial reference frame
%
% CONTRIBUTORS
%   Luca Deli
%
% VERSION
%   2026-01-03


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