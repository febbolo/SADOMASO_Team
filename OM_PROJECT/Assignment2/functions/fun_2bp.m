function dydt = fun_2bp(t, y, Mu_E)

%FUN_2BP Two-body problem equations of motion.
%
%   dydt = fun_2bp(t, y, Mu_E)
%
% PROTOTYPE
%   dydt = fun_2bp(t, y, Mu_E)
%
% DESCRIPTION
%   This function defines the equations of motion for the classical
%   two-body problem, describing the Keplerian motion of a spacecraft
%   orbiting a central body under the effect of its gravitational
%   attraction only.
%   The state vector is expressed in an inertial reference frame and
%   consists of position and velocity components. The acceleration is
%   computed assuming a point-mass central body.
%
% INPUT
%   t     [1x1] Time variable (not explicitly used, required by ODE solvers) [s]
%   y     [6x1] State vector:
%               y(1:3) = position vector                                   [km]
%               y(4:6) = velocity vector                                   [km/s]
%   Mu_E  [1x1] Gravitational parameter of the central body                [km^3/s^2]
%
% OUTPUT
%   dydt  [6x1] Time derivative of the state vector:
%               dydt(1:3) = velocity components                             [km/s]
%               dydt(4:6) = acceleration components                         [km/s^2]
%
% ASSUMPTIONS
%   - Two-body Keplerian motion (no perturbations).
%   - Point-mass central body.
%   - Inertial reference frame.
%
% CONTRIBUTORS
%   Luca Deli
%
% VERSION
%   2026-01-03


% -------------------- EQUATIONS OF MOTION --------------------
    r = y(1:3); % position vector
    v = y(4:6); % velocity vector
    
    % Gravitational acceleration
    a = -Mu_E * r / norm(r)^3;
    
    % State derivative
    dydt = [v; a];
end


