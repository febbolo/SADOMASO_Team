function [a,e,incl,raan,w,theta,h,Q] = car2kep(s,mu)
%CAR2KEP This functions returns the 6 classic orbital elements for a
% Keplerian orbit plus the semi-major axis and the direction cosines matrix

% INPUTS: s : state vector, written in the geocentric equatorial frame 
%         mu: planetary constant [km^3/s^2]

% OUTPUTS :  h : specific angular momentum [km^2/s]
%            e : eccentricity of the elliptical orbit [-]
%            incl : inclination of the orbit [deg]
%            raan : Right Ascension of the Ascending Node [deg]
%            w : argument of perigee    [deg]
%            theta : true anomaly [deg]
%            a : semi-major axis [km]
%            Q : direction cosines matrix (DCM) in general form : from
%            geocentric equatorial frame to perifocal frame. 
%
% AUTHORS : Amura Fabio


% Extract position and velocity vectors
r_v = s(1:3);
v_v = s(4:6);

r = norm(r_v);

v_r = dot(v_v,r_v)/r;

% Specific angular momentum
h_v = cross(r_v,v_v);
h = norm(h_v);

% Compute the inclination
incl = acosd(h_v(3)/h); 

% Avoid singularity: equatorial orbit
if(incl==0)
    N_hat = [1;0;0];
    raan = 0;
elseif(incl==180)
    N_hat = [-1;0;0];
    raan = 180;
else
    % Line of Nodes N
    N_v = cross([0,0,1],h_v);
    N = norm(N_v);
    N_hat = N_v/N;

    % Compute the Right Ascension of Ascending Node
    if(N_hat(2)>=0)
        raan = acosd(N_hat(1));
    else
        raan = 360 - acosd(N_hat(1));
    end
end

% Compute the eccentricity
e_v = 1/mu * (cross(v_v,h_v)) - r_v/r;
e = norm(e_v);

% Avoid singularity : circular orbit
if(e==0)
    e_hat = N_hat;
    w = 0;
else
    e_hat = e_v/e;
    if(e_v(3)>=0)
        w = acosd(dot(N_hat,e_hat));
    else
        w = 360 - acosd(dot(N_hat,e_hat));
    end
end

% Compute the true anomaly
if(v_r >= 0)
    theta = acosd(dot(e_hat, r_v / r));
else 
    theta = 360 - acosd(dot(e_hat, r_v / r));
end

% Compute the semi-major axis
a = h^2 / (mu * (1 - e^2));

% Compute the direction cosines matrix (DCM)
Q = [-sind(raan)*cosd(incl)*sind(w)+cosd(raan)*cosd(w),  cosd(raan)*cosd(incl)*sind(w)+sind(raan)*cosd(w),    sind(incl)*sind(w);...
     -sind(raan)*cosd(incl)*cosd(w)-cosd(raan)*sind(w),  cosd(raan)*cosd(incl)*cosd(w)-sind(raan)*sind(w),    sind(incl)*cosd(w);...
     sind(raan)*sind(incl),                           -cosd(raan)*sind(incl),                           cosd(incl)];

end