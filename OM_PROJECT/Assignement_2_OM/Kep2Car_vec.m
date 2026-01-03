function [r0, v0] = Kep2Car_vec(kep, mu)
%KEP2CAR_VEC Conversion from Keplerian elements to Cartesian state vector
%
%   [r0, v0] = Kep2Car_vec(kep, mu)
%
% PROTOTYPE
%   [r0, v0] = Kep2Car_vec(kep, mu)
%
% DESCRIPTION
%   This function converts a set of classical Keplerian orbital elements
%   into the corresponding Cartesian position and velocity vectors in an
%   inertial reference frame (ECI).
%   The transformation is performed by first expressing the state in the
%   perifocal reference frame (PQW) and then rotating it into the inertial
%   frame using the standard 3-1-3 Euler rotation sequence.
%
% INPUT
%   kep [6x1]  Vector of Keplerian elements:
%              kep(1) = a     Semi-major axis            [km]
%              kep(2) = e     Eccentricity               [-]
%              kep(3) = i     Inclination                [rad]
%              kep(4) = Omega Right ascension of AN      [rad]
%              kep(5) = omega Argument of pericenter     [rad]
%              kep(6) = theta True anomaly               [rad]
%
%   mu  [1x1]  Gravitational parameter of the central body [km^3/s^2]
%
% OUTPUT
%   r0 [3x1]   Position vector in inertial frame (ECI)     [km]
%   v0 [3x1]   Velocity vector in inertial frame (ECI)     [km/s]
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
%   2025-11-03
%

% -------------------------------------------------------------------------
% Unpack Keplerian elements
% -------------------------------------------------------------------------

a     = kep(1);    % semi-major axis
e     = kep(2);    % eccentricity
i     = kep(3);    % inclination
Raan  = kep(4);    % right ascension of ascending node
omega = kep(5);    % argument of pericenter
theta = kep(6);    % true anomaly

% -------------------------------------------------------------------------
% Orbital geometry
% -------------------------------------------------------------------------

p = a * (1 - e^2);         % semi-latus rectum
v = sqrt(mu / p);         % characteristic orbital velocity
r = p / (1 + e*cos(theta)); % radial distance

% -------------------------------------------------------------------------
% State expressed in perifocal (PQW) reference frame
% -------------------------------------------------------------------------

r_pqw = [ r*cos(theta);
          r*sin(theta);
          0 ];

v_pqw = v * [ -sin(theta);
               e + cos(theta);
               0 ];

% -------------------------------------------------------------------------
% Rotation matrix from PQW to inertial frame (ECI)
% 3-1-3 Euler rotation: Rz(Omega) * Rx(i) * Rz(omega)
% -------------------------------------------------------------------------

R = [ ...
 cos(Raan)*cos(omega) - sin(Raan)*sin(omega)*cos(i), ...
 -cos(Raan)*sin(omega) - sin(Raan)*cos(omega)*cos(i), ...
  sin(Raan)*sin(i);

 sin(Raan)*cos(omega) + cos(Raan)*sin(omega)*cos(i), ...
 -sin(Raan)*sin(omega) + cos(Raan)*cos(omega)*cos(i), ...
 -cos(Raan)*sin(i);

 sin(omega)*sin(i), ...
 cos(omega)*sin(i), ...
 cos(i) ];

% -------------------------------------------------------------------------
% Transformation to inertial Cartesian coordinates
% -------------------------------------------------------------------------

r0 = R * r_pqw;    % position vector in ECI
v0 = R * v_pqw;    % velocity vector in ECI

end