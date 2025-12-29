function T = T_RGT(k, m)

% T_RGT
%
% This function computes the orbital period required to obtain a repeating
% ground track for a spacecraft orbiting a rotating planet.
%
% The repeating ground track condition is defined by imposing that the
% spacecraft completes k orbital revolutions while the planet completes
% m rotations.
%
% INPUTS:
% k   Number of spacecraft orbital revolutions [-]
% m   Number of planet rotations during the same time interval [-]
%
% OUTPUT:
% T   Orbital period required for the repeating ground track [s]

% Planet rotation rate [rad/s]
We = deg2rad(15.04/3600);

% Orbital period for repeating ground track (radians are consistently used)
T = (2*pi*m)/(We*k);

end