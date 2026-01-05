function kep = car2kep(r, v, mu)
%CAR2KEP Conversion from Cartesian state to Keplerian elements.
%
%   kep = car2kep(r, v, mu)
%
% PROTOTYPE
%   kep = car2kep(r, v, mu)
%
% DESCRIPTION
%   This function converts a Cartesian state vector (position and velocity)
%   expressed in an inertial reference frame into the corresponding set of
%   classical Keplerian orbital elements.
%   The computation is based on the standard vector formulation involving
%   the specific angular momentum, node vector, and eccentricity vector.
%   Special cases such as circular and/or equatorial orbits are handled by
%   assigning conventional values to the undefined angular elements.
%
% INPUT
%   r   [3x1] Position vector in inertial reference frame                [km]
%   v   [3x1] Velocity vector in inertial reference frame                [km/s]
%   mu  [1x1] Gravitational parameter of the central body                [km^3/s^2]
%
% OUTPUT
%   kep [6x1] Vector of Keplerian elements:
%             kep(1) = a     Semi-major axis                             [km]
%             kep(2) = e     Eccentricity                                [-]
%             kep(3) = i     Inclination                                 [rad]
%             kep(4) = Omega Right ascension of the ascending node       [rad]
%             kep(5) = omega Argument of pericenter                      [rad]
%             kep(6) = f     True anomaly (or argument of latitude for
%                            circular orbits)                            [rad]
%
% ASSUMPTIONS
%   - Two-body Keplerian motion.
%   - Inertial reference frame centered at the central body.
%   - Elliptic orbits are expected (parabolic case handled but not nominal).
%   - For circular and/or equatorial orbits, undefined angular elements
%     are set according to standard conventions.
%
% CONTRIBUTORS
%   Luca Deli
%
% VERSION
%   2026-01-03

    r = r(:); 
    v = v(:);

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
