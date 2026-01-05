function [alpha, delta, lon_deg, lat_deg] = Ground_track(t, r, thetaG0, t0, wE)

%GROUND_TRACK Ground track computation from inertial position history.
%
%   [alpha, delta, lon_deg, lat_deg] = Ground_track(t, r, thetaG0, t0, wE)
%
% PROTOTYPE
%   [alpha, delta, lon_deg, lat_deg] = Ground_track(t, r, thetaG0, t0, wE)
%
% DESCRIPTION
%   This function computes the ground track (sub-satellite longitude and
%   latitude) of a spacecraft orbiting a rotating planet, starting from the
%   inertial position history r(t) expressed in an inertial reference frame
%   (e.g., ECI).
%   First, the right ascension (alpha) and declination (delta) of the
%   spacecraft are computed from the inertial position vector. Then, the
%   Greenwich sidereal angle is propagated from its initial value thetaG0
%   using the planet rotation rate wE, and the inertial angles are mapped to
%   body-fixed longitude and latitude. Longitudes are finally wrapped to
%   [-180, 180] degrees.
%
% INPUT
%   t       [Nx1] Time vector                                              [s]
%   r       [Nx3] Spacecraft position in inertial reference frame           [km]
%   thetaG0 [1x1] Greenwich sidereal angle at the initial time t0           [rad]
%   t0      [1x1] Initial time                                              [s]
%   wE      [1x1] Planet rotation rate                                      [rad/s]
%
% OUTPUT
%   alpha   [Nx1] Right ascension                                           [rad]
%   delta   [Nx1] Declination                                               [rad]
%   lon_deg [Nx1] Sub-satellite longitude (wrapped to [-180,180])           [deg]
%   lat_deg [Nx1] Sub-satellite latitude                                    [deg]
%
% ASSUMPTIONS
%   - r(t) is expressed in an inertial frame centered on the planet.
%   - Uniform planet rotation rate wE (constant angular velocity).
%   - The sub-satellite point is obtained neglecting Earth orientation
%     effects beyond a simple sidereal rotation (no polar motion, etc.).
%
% CONTRIBUTORS
%   Luca Deli
%
% VERSION
%   2026-01-03

N = length(t);
alpha = zeros(N,1);
delta = zeros(N,1);
lon = zeros(N,1);
lat = zeros(N,1);

for k = 1:N
    x = r(k,1);
    y = r(k,2);
    z = r(k,3);

    rnorm = norm([x y z]);

    % Declination
    delta(k) = asin(z / rnorm);

    % Right ascension
    alpha(k) = atan2(y, x);

    % Greenwich sidereal angle
    thetaG = thetaG0 + wE * (t(k) - t0);
    thetaG = mod(thetaG + pi, 2*pi) - pi;

    % Sub-satellite longitude and latitude
    lon(k) = alpha(k) - thetaG;
    lat(k) = delta(k);
end

% Convert to degrees and wrap longitude
lon_deg = rad2deg(lon);
lat_deg = rad2deg(lat);
lon_deg = wrapTo180(lon_deg);

end
