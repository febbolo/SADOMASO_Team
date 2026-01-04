function FLYBY = flyby_powered(vinf_min,vinf_plus,mu,R,hmin)
% FLYBY_POWERED This function solves the hyperbola in fly-by flight dynamics
% application by exploiting geometrical formulas of the hyperbolas uniquely
% defined by vinf, Delta.
%
% INPUT : vinf_min : Relative Velocity of S/C wrt planet vector before powered fly-by
%                 in Heliocentric Ecliptic Inertial Frame [km/s]
%         vinf_plus : Relative Velocity of S/C wrt planet vector after powered fly-by
%                 in Heliocentric Ecliptic Inertial Frame [km/s]
%         mu : Planetary constant [km^3/s^2]
%         R : Radius of the planet [km^3/s^2]
%         hmin : Altitude of the atmosphere (if present) of the planet
%
% OUTPUT : FLYBY: Struct containing all the data for the two hyperbolas
%          FLYBY Struct Description:
%          FLYBY.HYP1   : Struct containing parameters for the first hyperbola
%          .vinf        : Relative velocity of the spacecraft before the flyby [km/s]
%          .a_hyp       : Semi-major axis of the first hyperbola [km]
%          .e_hyp       : Eccentricity of the first hyperbola
%          .delt        : Delta angle for the first hyperbola [degrees]
%          .Delta1      : Impact parameter for the first hyperbola
%           .vp_hyp     : Periapsis velocity of the first hyperbola [km/s]
%           .vp_hyp_v   : Periapsis velocity vector of the first hyperbola
%           .h_hyp      : Angular momentum of the first hyperbola
%           .h_hyp_v    : Angular momentum vector of the first hyperbola
%           .theta_inf  : Orientation angle of the first hyperbola [degrees]
%
%           FLYBY.HYP2  : Struct containing parameters for the second hyperbola
%           .vinf       : Relative velocity of the spacecraft after the flyby [km/s]
%           .a_hyp      : Semi-major axis of the second hyperbola [km]
%           .e_hyp      : Eccentricity of the second hyperbola
%           .delta2     : Delta angle for the second hyperbola [degrees]
%           .Delta2     : Impact parameter for the second hyperbola
%           .vp_hyp     : Periapsis velocity of the second hyperbola [km/s]
%           .vp_hyp_v   : Periapsis velocity vector of the second hyperbola
%           .h_hyp      : Angular momentum of the second hyperbola
%           .h_hyp_v    : Angular momentum vector of the second hyperbola
%           .theta_inf  : Orientation angle of the second hyperbola [degrees]
%
%           FLYBY.delta       : Angle between the two velocity vectors [degrees]
%           FLYBY.rp          : Radius of periapsis [km]
%           FLYBY.rp_v        : Radius of periapsis vector
%           FLYBY.Deltavp     : Change in periapsis velocity [km/s]
%           FLYBY.Deltavp_v   : Change in periapsis velocity vector
%           FLYBY.u_hat       : Unit vector normal to the plane of the hyperbolas
%           FLYBY.Deltav      : Change in velocity of the spacecraft [km/s]
%           FLYBY.vinf_min    : Initial relative velocity vector [km/s]
%           FLYBY.vinf_plus   : Final relative velocity vector [km/s]
%           FLYBY.flagtype    : Type of flyby ("Powered" or "Unpowered")
%
% 
% AUTHORS :     Amura Fabio
%
%


% Find vinf 1
vinf1 = norm(vinf_min); %[km/s]
% Find vinf 2
vinf2 = norm(vinf_plus); %[km/s]


% Find semi-major axis of hyperbola 1
a_hyp1 = -mu / (vinf1)^2; %[km]
% Find semi-major axis of hyperbola 2
a_hyp2 = -mu / (vinf2)^2; %[km]

% Compute delta
delta = acosd(dot(vinf_min,vinf_plus) / (norm(vinf_min)*norm(vinf_plus)));

% Compute rp
rp = find_rp(vinf1,vinf2,delta,mu,R,hmin);

% Find eccentricity of hyperbola 1
e_hyp1 = 1 + (rp*(vinf1)^2)/mu;
% Find eccentricity of hyperbola 2
e_hyp2 = 1 + (rp*(vinf2)^2)/mu;

% Compute delta of hyperbola 1
delta1 = 2 * asind(1/e_hyp1); %[deg]
% Compute delta of hyperbola 2
delta2 = 2 * asind(1/e_hyp2); %[deg]

% Compute vp of hyperbola 1
vp_hyp1 = sqrt( (vinf1)^2 + (2*mu/rp) ); %[km/s]
% Compute vp of hyperbola 2
vp_hyp2 = sqrt( (vinf2)^2 + (2*mu/rp) ); %[km/s]

% Find plane of the hyperbola and angular momentum vector
u_hat = cross(vinf_min,vinf_plus) / norm(cross(vinf_min,vinf_plus));
h_hyp1 = sqrt(mu * a_hyp1 * (1-e_hyp1^2));
h_hyp_v1 = h_hyp1 * u_hat;
h_hyp2 = sqrt(mu * a_hyp2 * (1-e_hyp2^2));
h_hyp_v2 = h_hyp2 * u_hat;

%Find orientation of the apse line using Rodrigues' formula
theta_inf1 = acosd(-1/e_hyp1);
theta_inf2 = acosd(-1/e_hyp2);
beta1 = 180-theta_inf1;
e_hyp_hat = vinf_min*cosd(-beta1) + cross(u_hat,vinf_min) * sind(-beta1) ...
            + u_hat * dot(u_hat,vinf_min) * (1 - cosd(-beta1));
e_hyp_hat = e_hyp_hat / norm(e_hyp_hat);

% Incoming hyperbola plane
t_hat_in = vinf_min / norm(vinf_min);
h_hat_in = cross(e_hyp_hat, t_hat_in);
h_hat_in = h_hat_in / norm(h_hat_in);
p_hat_in = cross(h_hat_in, e_hyp_hat);
p_hat_in = p_hat_in / norm(p_hat_in);

% Outgoing hyperbola plane/basis 
t_hat_out = vinf_plus / norm(vinf_plus);
h_hat_out = cross(e_hyp_hat, t_hat_out);
h_hat_out = h_hat_out / norm(h_hat_out);
p_hat_out = cross(h_hat_out, e_hyp_hat);
p_hat_out = p_hat_out / norm(p_hat_out);

% Periapsis velocities 
vp_hyp1_v = vp_hyp1 * p_hat_in;
vp_hyp2_v = vp_hyp2 * p_hat_out;

% Find Deltavp_v
Deltavp_v = vp_hyp2_v - vp_hyp1_v;
Deltavp   = norm(Deltavp_v);

% Find Deltav
Deltav_v = vinf_plus - vinf_min;
Deltav = norm(Deltav_v);

% Compute Delta (Impact parameter) of hyperbola 1
Delta1 = -a_hyp1*e_hyp1*cosd(delta1/2);
% Compute Delta (Impact parameter) of hyperbola 2
Delta2 = -a_hyp2*e_hyp2*cosd(delta2/2);

% Computing radius of pericentre vector
rp_hyp = rp * e_hyp_hat;

% Constructing struct FLYBY : hyperbola 1
FLYBY.HYP1.vinf = vinf1;
FLYBY.HYP1.a_hyp = a_hyp1;
FLYBY.HYP1.e_hyp = e_hyp1;
FLYBY.HYP1.delta1 = delta1;
FLYBY.HYP1.Delta1 = Delta1;
FLYBY.HYP1.vp_hyp = vp_hyp1;
FLYBY.HYP1.vp_hyp_v = vp_hyp1_v;
FLYBY.HYP1.h_hyp = h_hyp1;
FLYBY.HYP1.h_hyp_v = h_hyp_v1;
FLYBY.HYP1.theta_inf = theta_inf1;
% Constructing struct FLYBY : hyperbola 2
FLYBY.HYP2.vinf = vinf2;
FLYBY.HYP2.a_hyp = a_hyp2;
FLYBY.HYP2.e_hyp = e_hyp2;
FLYBY.HYP2.delta2 = delta2;
FLYBY.HYP2.Delta2 = Delta2;
FLYBY.HYP2.vp_hyp = vp_hyp2;
FLYBY.HYP2.vp_hyp_v = vp_hyp2_v;
FLYBY.HYP2.h_hyp = h_hyp2;
FLYBY.HYP2.h_hyp_v = h_hyp_v2;
FLYBY.HYP2.theta_inf = theta_inf2;
% Constructing struct FLYBY : general
FLYBY.delta = delta;
FLYBY.rp = rp;
FLYBY.rp_v = rp_hyp;
FLYBY.Deltavp = Deltavp;
FLYBY.Deltavp_v = Deltavp_v;
FLYBY.u_hat = u_hat; 
FLYBY.Deltav_flyby = Deltav;
FLYBY.vinf_min = vinf_min;
FLYBY.vinf_plus = vinf_plus;

% Flag : powered or unpowered flyby
if(vinf1 ~= vinf2)
    FLYBY.flagtype = "Powered";
else
    FLYBY.flagtype = "Unpowered";
end


end



