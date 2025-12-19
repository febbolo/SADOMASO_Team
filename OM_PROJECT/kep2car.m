function [s,Q] = kep2car(a,e,incl,raan,w,theta,mu)

%KEP2CAR This functions returns the state vector in the geocentric equatorial frame 
% given in input the 6 classical orbital elements for a
% Keplerian orbit plus the planetary constant

% INPUTS:    a : semi-major axis of the orbit [km]
%            e : eccentricity of the elliptical orbit [-]
%            incl : inclination of the orbit [deg]
%            raan : Right Ascension of the Ascending Node [deg]
%            w : argument of perigee    [deg]
%            theta : true anomaly [deg]
%            mu: planetary constant [km^3/s^2]

% OUTPUTS :  s : state vector, written in the geocentric equatorial frame 
%            Q : direction cosines matrix (DCM) from perifocal frame to geocentric equatorial frame  

% Computing h
h = sqrt(mu * a * (1 - e^2));

r_v_perifocal = h^2/(mu*(1+e*cosd(theta))) * [cosd(theta);...
                                              sind(theta);...
                                              0];
v_v_perifocal = mu/h * [-sind(theta);...
                        e+cosd(theta);...
                        0];

%DCM from geocentric to perifocal
Q = [-sind(raan)*cosd(incl)*sind(w)+cosd(raan)*cosd(w),  cosd(raan)*cosd(incl)*sind(w)+sind(raan)*cosd(w),    sind(incl)*sind(w);...
     -sind(raan)*cosd(incl)*cosd(w)-cosd(raan)*sind(w),  cosd(raan)*cosd(incl)*cosd(w)-sind(raan)*sind(w),    sind(incl)*cosd(w);...
     sind(raan)*sind(incl),                           -cosd(raan)*sind(incl),                           cosd(incl)];

%Q'Q = I -> DCM orthogonal, so take Q' to go from perifocal to geocentric
%equatorial
Q_T = Q';
s = [Q_T * r_v_perifocal; Q_T * v_v_perifocal];

end