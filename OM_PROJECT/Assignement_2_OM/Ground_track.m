function [alpha, delta, lon_deg, lat_deg] = Ground_track(t, r, thetaG0, t0, wE)
% GROUND_TRACK
%
% This function computes the ground track of a spacecraft orbiting a planet.
% Starting from the inertial position history r(t), the right ascension and
% declination are computed and then converted into longitude and latitude on
% the rotating planet through the Greenwich sidereal angle.
%
% INPUTS:
% t        Time vector [s]
% r        Spacecraft position in the inertial reference frame [km]
% thetaG0  Greenwich sidereal angle at the initial time t0 [rad]
% t0       Initial time [s]
%
% OUTPUTS:
% alpha    Right ascension of the spacecraft [rad]
% delta    Declination of the spacecraft [rad]
% lon_deg      Longitude of the sub-satellite point [deg]
% lat_deg      Latitude of the sub-satellite point [deg]
% Rotation rate [rad/s]

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

% Converting in degrees and wrapping
lon_deg = rad2deg(lon);
lat_deg = rad2deg(lat);
lon_deg = wrapTo180(lon_deg);

end