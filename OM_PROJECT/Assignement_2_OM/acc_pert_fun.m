function acc_pert_vec = acc_pert_fun(t, s, parameters)
%ACC_PERT_FUN Perturbing acceleration model (J2 + Moon third-body) in ECI.
%
%   acc_pert_vec = acc_pert_fun(t, s, parameters)
%
% PROTOTYPE
%   acc_pert_vec = acc_pert_fun(t, s, parameters)
%
% DESCRIPTION
%   This function computes the total perturbing acceleration acting on a
%   spacecraft in Earth-centered inertial coordinates (ECI). The model
%   includes:
%   (i) the acceleration due to Earth oblateness (J2), and
%   (ii) the Moon third-body gravitational perturbation, computed using the
%   differential acceleration formulation (spacecraft-to-Moon minus
%   Earth-to-Moon).
%   The Moon position is obtained from ephemerides through ephMoon at the
%   current epoch, updated from an initial MJD2000 time by adding the
%   elapsed integration time.
%
% INPUT
%   t           [1x1] Time from initial epoch                               [s]
%   s           [6x1] Cartesian state vector in ECI:
%                     s(1:3) = position vector                              [km]
%                     s(4:6) = velocity vector                              [km/s]
%   parameters  [1x5] Model parameters vector:
%                     parameters(1) = J2                                    [-]
%                     parameters(2) = mu (Earth)                             [km^3/s^2]
%                     parameters(3) = mass (Earth)                           [kg]  (not used)
%                     parameters(4) = Radius (Earth)                         [km]
%                     parameters(5) = t0_mjd2000 initial epoch               [days]
%
% OUTPUT
%   acc_pert_vec [3x1] Total perturbing acceleration in ECI                  [km/s^2]
%
% ASSUMPTIONS
%   - ECI inertial frame centered at Earth.
%   - Earth gravity field modelled only through the J2 zonal harmonic.
%   - The Moon is treated as a point mass; its position is obtained from
%     ephemerides (ephMoon) referred to the same ECI frame.
%   - The third-body acceleration is computed as a differential term to
%     remove the common Earth-Moon acceleration (indirect term).
%
% CONTRIBUTORS
%   Luca Deli
%
% VERSION
%   2026-01-03

% parameters in order:
% parameters = [Primary.J2, Primary.mu, Primary.mass, Primary.Radius, t0_mjd2000]

Primary.J2     = parameters(1);
Primary.mu     = parameters(2);
Primary.mass   = parameters(3); %#ok<NASGU> % not used in this function
Primary.Radius = parameters(4);
t0_mjd2000     = parameters(5);

r = s(1:3);

x = r(1);
y = r(2);
z = r(3);

rnorm = norm(r);

% ____________________ 3 BODY PERTURBATIONS (Moon) ___________________

% Update epoch in MJD2000 (seconds -> days)
t_mjd2000 = t0_mjd2000 + t/86400;

moon.mu = astroConstants(20);

[x3Body, ~] = ephMoon(t_mjd2000); % Moon position in ECI [km]

RSC_3B      = x3Body' - r;          % spacecraft -> Moon
RSC_3B_norm = norm(RSC_3B);

RCB_3B      = x3Body';              % Earth -> Moon
RCB_3B_norm = norm(RCB_3B);

% Acceleration due to third body (Moon), differential formulation
a_3bp_moon = moon.mu * ( (RSC_3B/(RSC_3B_norm^3)) - (RCB_3B/(RCB_3B_norm^3)) );

% _______________________ J2 PERTURBATION __________________________

a_J2 = ((3/2) * Primary.mu * (Primary.Radius^2) * Primary.J2 / (rnorm^4)) * ...
       [ (x/rnorm) * (5*(z^2/rnorm^2) - 1);
         (y/rnorm) * (5*(z^2/rnorm^2) - 1);
         (z/rnorm) * (5*(z^2/rnorm^2) - 3) ];

% Total perturbing acceleration in ECI
acc_pert_vec = a_J2 + a_3bp_moon;

end
