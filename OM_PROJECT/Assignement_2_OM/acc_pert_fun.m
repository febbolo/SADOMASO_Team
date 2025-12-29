function acc_pert_vec = acc_pert_fun(t, s, parameters)

% parameters in order:
% parameters = [Primary.J2, Primary.mu, Primary.mass, Primary.Radius, t0_mjd2000]

%           function
%           acc_pert_fun is used to compute the perturbing accelerations.
%           The perturbations considered are J2 and the Moon third-body effect.
%
% PROTOTYPE
%           acc_pert_vec = acc_pert_fun(t, s, parameters)
%
% INPUT
%           t [s]        time from the beginning of iterations
%           s [6x1]      state vector [r; v] in ECI
%           parameters   constants and initial epoch (MJD2000)
%
% OUTPUT
%           acc_pert_vec [3x1] total perturbing acceleration in ECI [km/s^2]

Primary.J2 = parameters(1);
Primary.mu = parameters(2);
Primary.mass = parameters(3);
Primary.Radius = parameters(4);
t0_mjd2000 = parameters(5);

r = s(1:3);

x = r(1);
y = r(2);
z = r(3);

rnorm = norm(r);

% ____________________ 3 BODY PERTURBATIONS (Moon) ___________________

% Every iteration I need to update the julian date adding seconds in day
% format in a way to find each ttime the right Ephemerides with the function

t_mjd2000 = t0_mjd2000 + t/86400; % 60*60*24 is a day in seconds

moon.mu = astroConstants(20);

[x3Body, ~] = ephMoon(t_mjd2000);     % Moon position in ECI (km)

RSC_3B = x3Body' - r;  % vectorial distance between spacecraft and third body
RSC_3B_norm = norm(RSC_3B);  % distance in norm

RCB_3B = x3Body'; % vectorial distance between third body and earth: I find it throught ephemerides
RCB_3B_norm = norm(RCB_3B);

% Acceleration due to third body (Moon)
a_3bp_moon = moon.mu * ( (RSC_3B/(RSC_3B_norm^3)) - (RCB_3B/(RCB_3B_norm^3)) );

% _______________________ J2 PERTURBATION __________________________
% Calculate the total perturbing acceleration due to J2 in vecctorial form
a_J2 = ((3/2) * Primary.mu * (Primary.Radius^2) * Primary.J2 / (rnorm^4)) * ...
       [ (x/rnorm) * (5*(z^2/rnorm^2) - 1);
         (y/rnorm) * (5*(z^2/rnorm^2) - 1);
         (z/rnorm) * (5*(z^2/rnorm^2) - 3) ];

% Total perturbing acceleration in ECI
acc_pert_vec = a_J2 + a_3bp_moon;

end
