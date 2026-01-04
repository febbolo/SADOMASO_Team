function rp = find_rp(vinf_min, vinf_plus, delta, mu, R, h_min)
% FIND_RP  This function computes the pericenter radius rp that produces a turn angle equal to delta.
%  See also FLYBY_POWERED function
% 
% INPUT :   vinf_minus : incoming hyperbolic excess speed [km/s]
%           vinf_plus : outgoing hyperbolic excess speed [km/s]
%           mu: gravitational parameter of the planet [km^3/s^2]
%           delta : turn angle, delta = delta1 + delta2 (powered FB) [deg]
%           R : radius of the planet [km]
%           hmin : Minimum altitude due to atmosphere of planet [km]
%
% OUTPUT : rp : radius of pericentre of both hyperbolic legs [km]
% 
%
% AUTHORS :     Amura Fabio
%

% Computing radius of the planet

% Eccentricities as functions of rp
e_min = @(rp) 1 + (rp .* vinf_min.^2) / mu;
e_plus  = @(rp) 1 + (rp .* vinf_plus.^2) / mu;

% Partial deflection angles
delta_min = @(rp) 2 .* asind(1 ./ e_min(rp)); % [deg]
delta_plus  = @(rp) 2 .* asind(1 ./ e_plus(rp));  % [deg]

% Total deflection angle (average of the two)
delta_total = @(rp) (delta_min(rp) + delta_plus(rp)) / 2;

% Function to solve: delta_total(rp) - delta = 0
f = @(rp) delta_total(rp) - delta;

% Solve using fzero : automatically neglecting rp < R+hmin
rpmin = R + h_min; % [km]
rpmax = 1e12; %[km]

% Check feasability 
if f(rpmin) * f(rpmax) > 0
    rp = NaN;        % IMPOSSIBLE FLYBY
    return
end

rp = fzero(f, [rpmin rpmax]);

end