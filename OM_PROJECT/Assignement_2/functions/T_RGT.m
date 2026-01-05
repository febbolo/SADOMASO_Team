function T = T_RGT(k, m)
%T_RGT Orbital period for repeating ground track condition.
%
%   T = T_RGT(k, m)
%
% PROTOTYPE
%   T = T_RGT(k, m)
%
% DESCRIPTION
%   This function computes the orbital period required to achieve a
%   repeating ground track for a spacecraft orbiting a uniformly rotating
%   planet. The repeating ground track condition is enforced by requiring
%   that the spacecraft completes k orbital revolutions while the planet
%   completes m full rotations over the same time interval.
%   Under this assumption, the required orbital period is obtained from the
%   commensurability between the mean orbital motion and the planet rotation
%   rate.
%
% INPUT
%   k [1x1] Number of spacecraft orbital revolutions                    [-]
%   m [1x1] Number of planet rotations in the same time interval        [-]
%
% OUTPUT
%   T [1x1] Orbital period required for repeating ground track           [s]
%
% ASSUMPTIONS
%   - Uniform planet rotation rate (constant angular velocity).
%   - Two-body mean motion representation (no precession effects included).
%   - The repeating condition is defined only by the k:m commensurability.
%
% CONTRIBUTORS
%   Luca Deli
%
% VERSION
%   2026-01-03

% Planet rotation rate [rad/s]
We = deg2rad(15.04/3600);

% Orbital period for repeating ground track
T = (2*pi*m)/(We*k);

end
