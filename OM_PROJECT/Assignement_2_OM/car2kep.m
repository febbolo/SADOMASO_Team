function kep = car2kep(r, v, mu)
%CAR2KEP Convert Cartesian state to classical Keplerian elements.
%
% INPUT:
%   r  [3x1] position (km)
%   v  [3x1] velocity (km/s)
%   mu [1x1] gravitational parameter (km^3/s^2)
%
% OUTPUT:
%   kep = [a; e; i; Omega; omega; f]
%     a     semi-major axis (km)
%     e     eccentricity (-)
%     i     inclination (rad)
%     Omega RAAN (rad)
%     omega argument of perigee (rad)
%     f     true anomaly (rad)

    r = r(:); v = v(:);

    rnorm = norm(r);
    vnorm = norm(v);

    hvec = cross(r, v);
    h    = norm(hvec);

    k = [0;0;1];
    nvec = cross(k, hvec);
    n    = norm(nvec);

    evec = (1/mu) * ( cross(v, hvec) - mu * r / rnorm );
    e    = norm(evec);

    % specific mechanical energy
    eps = vnorm^2/2 - mu/rnorm;

    % semi-major axis
    if abs(e-1) > 1e-10
        a = -mu/(2*eps);
    else
        a = Inf; % parabolic (not expected here)
    end

    % inclination
    i = acos( hvec(3)/h );

    % RAAN Omega
    if n > 1e-12
        Omega = atan2(nvec(2), nvec(1));
    else
        Omega = 0; % equatorial
    end

    % argument of perigee omega
    if (n > 1e-12) && (e > 1e-12)
        omega = atan2( dot(cross(nvec, evec), hvec)/h, dot(nvec, evec) );
    else
        omega = 0; % circular or equatorial: undefined
    end

    % true anomaly f
    if e > 1e-12
        f = atan2( dot(cross(evec, r), hvec)/h, dot(evec, r) );
    else
        % circular: use argument of latitude u instead
        if n > 1e-12
            u = atan2( dot(cross(nvec, r), hvec)/h, dot(nvec, r) );
        else
            u = atan2(r(2), r(1)); % circular equatorial
        end
        f = u; % in circular case f not defined; we return u
    end

    % normalize angles to [0, 2pi)
    Omega = mod(Omega, 2*pi);
    omega = mod(omega, 2*pi);
    f     = mod(f, 2*pi);

    kep = [a; e; i; Omega; omega; f];
end