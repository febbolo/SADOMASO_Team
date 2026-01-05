function dkep = eq_motion_Gauss_J2_Moon(t, kep, acc_pert_fun, parameters)

%EQ_MOTION_GAUSS_J2_MOON Gauss planetary equations (RSW) for perturbed orbits.
%
%   dkep = eq_motion_Gauss_J2_Moon(t, kep, acc_pert_fun, parameters)
%
% PROTOTYPE
%   dkep = eq_motion_Gauss_J2_Moon(t, kep, acc_pert_fun, parameters)
%
% DESCRIPTION
%   This function defines the Gauss planetary equations for the propagation
%   of classical Keplerian elements under generic perturbations.
%   The perturbing acceleration is provided in the inertial frame (ECI)
%   through the function handle acc_pert_fun and is projected into the local
%   Radial-Transverse-Out-of-plane (RSW) frame built from the instantaneous
%   osculating state (r,v). The resulting RSW components (a_r, a_s, a_w) are
%   used to compute the time derivatives of the Keplerian elements.
%   The function is intended for perturbation models including, e.g., Earth
%   oblateness (J2) and third-body gravity (Moon/Sun), as implemented inside
%   acc_pert_fun.
%
% INPUT
%   t           [1x1] Time from initial epoch                               [s]
%   kep         [6x1] Keplerian elements:
%                     kep(1) = a     Semi-major axis                        [km]
%                     kep(2) = e     Eccentricity                           [-]
%                     kep(3) = i     Inclination                            [rad]
%                     kep(4) = Raan  Right ascension of AN                  [rad]
%                     kep(5) = omega Argument of pericenter                 [rad]
%                     kep(6) = f     True anomaly                           [rad]
%
%   acc_pert_fun [func] Function handle providing the perturbing acceleration:
%                     acc_ECI = acc_pert_fun(t, s_cart, parameters)
%                     with s_cart = [r; v] in ECI                           [km; km/s]
%
%   parameters  [1xN] Vector of model parameters (user-defined), e.g.:
%                     [J2, mu, mass, Radius, t0_mjd2000]
%
% OUTPUT
%   dkep        [6x1] Time derivatives of Keplerian elements:
%                     [da/dt; de/dt; di/dt; dRaan/dt; domega/dt; df/dt]
%
% ASSUMPTIONS
%   - Osculating Keplerian elements (instantaneous orbit).
%   - Perturbations are treated as small accelerations added to the two-body dynamics.
%   - The RSW frame is built from the instantaneous state (r,v).
%   - Elliptic orbit (0 <= e < 1). Note: omega_dot and f_dot include 1/e terms.
%   - inc ~= 0 to avoid singularity in dRaan/dt and omega_dot (1/sin(inc)).
%
% CONTRIBUTORS
%   Luca Deli
%
% VERSION
%   2026-01-03

a     = kep(1);
e     = kep(2);
inc   = kep(3);
Raan  = kep(4);
omega = kep(5);
f     = kep(6);

Primary.mu = parameters(2);

% Convert kep -> cartesian to build the local RSW frame
[r, v] = Kep2Car_vec(kep, Primary.mu);
s = [r; v];

% Perturbing acceleration in ECI
acc_ECI = acc_pert_fun(t, s, parameters);

% Orbital relations
p     = a * (1 - e^2);
rnorm = p / (1 + e*cos(f));
hvec  = cross(r, v);
h     = norm(hvec);

% RSW frame unit vectors
r_rsw = r / rnorm;
w_rsw = hvec / h;
s_rsw = cross(w_rsw, r_rsw);

% Projecting the accelerations into RSW
a_r = dot(acc_ECI, r_rsw);
a_s = dot(acc_ECI, s_rsw);
a_w = dot(acc_ECI, w_rsw);

% --- Gauss planetary equations (RSW form) ---
a_dot   = (2*a^2/h) * ( e*sin(f)*a_r + (p/rnorm)*a_s );
e_dot   = (1/h) * ( p*sin(f)*a_r + ((p + rnorm)*cos(f) + rnorm*e)*a_s );
inc_dot = (rnorm*cos(f + omega)/h) * a_w;

Raan_dot  = (rnorm*sin(f + omega)/(h*sin(inc))) * a_w;
omega_dot = (1/(h*e)) * ( -p*cos(f)*a_r + (p + rnorm)*sin(f)*a_s ) ...
          - (rnorm*sin(f + omega)*cos(inc)/(h*sin(inc))) * a_w;

f_dot   = (h/rnorm^2) + (1/(h*e)) * ( p*cos(f)*a_r - (p + rnorm)*sin(f)*a_s );

% --- Pack output ---
dkep = [a_dot; e_dot; inc_dot; Raan_dot; omega_dot; f_dot];

end
