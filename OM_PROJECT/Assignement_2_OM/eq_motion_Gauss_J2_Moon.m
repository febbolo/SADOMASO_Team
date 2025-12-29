function dkep = eq_motion_Gauss_J2_Moon(t, kep, acc_pert_fun, parameters)

%           function
%           Gauss planetary equations (RSW) for perturbed orbits.
%           Perturbations are provided by acc_pert_fun in ECI and projected in RSW.
%
% PROTOTYPE
%           dkep = eq_motion_Gauss_J2_Moon(t, kep, acc_pert_fun, parameters)
%
% INPUT
%           t [s]                 time from initial epoch
%           kep [6x1]             [a; e; i; Raan; omega; f] (a in km, angles in rad)
%           acc_pert_fun handle   a_pert_ECI = acc_pert_fun(t, s_cart, parameters)
%           parameters            [J2, mu, mass, Radius, t0_mjd2000]
%
% OUTPUT
%           dkep [6x1]            time derivatives of Keplerian elements

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

    % --- Gauss planetary equations (RSW form), I define derivates of orbital parameters ---

    % derivative of semi-major axis in time
    a_dot     = (2*a^2/h) * ( e*sin(f)*a_r + (p/rnorm)*a_s );

    % derivative of eccentricity in time
    e_dot     = (1/h) * ( p*sin(f)*a_r + ((p + rnorm)*cos(f) + rnorm*e)*a_s );

    % derivative of inclination in time
    inc_dot   = (rnorm*cos(f + omega)/h) * a_w;

    % derivative of anomaly of pericentre in time
    Raan_dot = (rnorm*sin(f + omega)/(h*sin(inc))) * a_w;

    % derivative of Raan in time
    omega_dot = (1/(h*e)) * ( -p*cos(f)*a_r + (p + rnorm)*sin(f)*a_s ) - (rnorm*sin(f + omega)*cos(inc)/(h*sin(inc))) * a_w;

    % derivative of f0 in time
    f_dot     = (h/rnorm^2) + (1/(h*e)) * ( p*cos(f)*a_r - (p + rnorm)*sin(f)*a_s );

    % --- Pack output ---
    dkep = [a_dot; e_dot; inc_dot; Raan_dot; omega_dot; f_dot];

end